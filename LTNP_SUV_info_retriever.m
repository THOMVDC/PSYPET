

%tracer='MK6240';
petdir = '/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/MK62/output/psypet_cat12_output/*0*';
rawdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/MK62_4s';
outputexcel = '/Volumes/LaCie/Thomas/Projects/MK_paper/SUV_info.xlsx';

l = dir(petdir);
T = cell(length(l),6);
T{1,1} = 'subject';
T{1,1+1} =  'bodyWeight';
T{1,1+2} =  'injdosis';
T{1,1+3} =  'acqtime';
T{1,1+4} =  'injtime';
T{1,1+5} =  'halftime';
T{1,1+6} =  'rescale_factor';

for i = 1:numel(l)
    
    % define subject names based on dicom directory and convert to SUV
    subject = l(i).name;
    %psypet=dir(fullfile(l(i).folder,subject,'psypet*'));
    %psypet=psypet(1).name;
    %dicomdir = fullfile(l(i).folder,subject,psypet,tracer,'DCM');
    dicomdir=fullfile(rawdir,subject);
    
    % outname=['SUV_' subject];
    %[ptname,ptID,bodyweight, injdosis, acqtime, injtime, halftime]=LTNP_dcm2SUV(dicomdir,outputdir,outname)
    
    % Loading dicom header
    list=dir(dicomdir);
    list=list(~startsWith({list.name}, '.')); % removing hidden files from the list
    dicomimg=fullfile(dicomdir,list(10).name); % grabs tenth image (pseudorandom)
    info = dicominfo(dicomimg,'UseDictionaryVR',true);
    
    if isfield(info.PatientName,'FamilyName') && isfield(info.PatientName,'GivenName') 
        ptname=[info.PatientName.FamilyName '_' info.PatientName.GivenName];
    else
        ptname='unknown, please refer to ptID';
    end
    
    if isfield(info,'PatientID')
        ptID=info.PatientID;
    else
        ptID='unknown, please refer to ptname';
    end
    
    if isfield(info,'PatientWeight')
        bodyweight = 1000*info.PatientWeight; % in g
    else
        error('no PatientWeight specified in dcm')
    end
    
    
    injdosis=0;
    while injdosis==0
        if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose')
            if info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose ~= 0
                injdosis = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % in Bq
            else
                info.RadiopharmaceuticalInformationSequence.Item_1 = rmfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose');
            end
        elseif isfield(info,'Private_0009_1038')
            if info.Private_0009_1038 ~= 0
                injdosis = 1000000*info.Private_0009_1038;
            else
                info=rmfield(info,'Private_0009_1038');
            end
        elseif isfield(info,'Private_0009_103a')
            if info.Private_0009_103a ~= 0
                injdosis = 1000000*info.Private_0009_103a;
            else
                info=rmfield(info,'Private_0009_103a');
            end
        else
            error('no RadionuclideTotalDose specified in dcm')
        end
    end
    
%     if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose')
%         injdosis = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % in Bq
%         display='Radiopharm';
%         if injdosis <= 0
%             error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
%         end
%     elseif isfield(info,'Private_0009_1038')
%         injdosis = 1000000*info.Private_0009_1038;
%         display='Private';
%         if injdosis <= 0
%             error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
%         end
%     else
%         error('no RadionuclideTotalDose specified in dcm')
%     end
       
    if isfield(info,'AcquisitionTime')
        acqtime= info.AcquisitionTime;
        if length(acqtime)>6
            acqtime=acqtime(1:6);
        end
    else
        error('no AcquisitionTime specified in dcm')
    end
    
    if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadiopharmaceuticalStartTime')
        injtime= info.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
        if length(injtime)==6
            injtime=[injtime '.000'];
        end
    else
        error('no RadiopharmaceuticalStartTime specified in dcm')
    end
    
    if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideHalfLife')
        halftime=info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife; % in seconds
    else
        error('no RadionuclideHalfLife specified in dcm')
    end
    
    % Calculate rescale factor
    dt=datetime(acqtime,'InputFormat','HHmmss')-datetime(injtime,'InputFormat','HHmmss.SSS');
    dt=seconds(dt); % in seconds
    if dt<3000
        error('no RadionuclideHalfLife specified in dcm')
    end
    
    rescale_factor=(bodyweight/injdosis) * (1/2)^(-dt/halftime);
    
    % write dcm info into tabel to check missing data
    T{1+i,1} = subject;
    T{1+i,1+1} = bodyweight;
    T{1+i,1+2} = injdosis;
    T{1+i,1+3} = acqtime;
    T{1+i,1+4} = injtime;
    T{1+i,1+5} = halftime;
    T{1+i,1+6} = rescale_factor;
    
end

writecell(T,outputexcel);


% Convert dcm2niix
% requires: pet dcm in subject folders (one deep?)
% next step: extract brain based on CT 

tracer='UCBJ';
petdir = '/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/UCBJ_4s';
outputdir = '/Users/tvdcas1/Desktop/';

l = dir(petdir);
table = cell(length(l),6);
table{1,1} = 'subject';
table{1,1+1} =  'bodyWeight';
table{1,1+2} =  'injdosis';
table{1,1+3} =  'acqtime';
table{1,1+4} =  'injtime';
table{1,1+5} =  'halftime';

for i = 51:numel(l)
    
    % define subject names based on dicom directory and convert to SUV
    subject = l(i).name;
    psypet=dir(fullfile(l(i).folder,subject,'psypet*'));
    psypet=psypet(1).name;
    dicomdir = fullfile(l(i).folder,subject,psypet,tracer,'DCM');
    
    % outname=['SUV_' subject];
    %[ptname,ptID,bodyweight, injdosis, acqtime, injtime, halftime]=LTNP_dcm2SUV(dicomdir,outputdir,outname)

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
    
    if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose')
        injdosis = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % in Bq
        display='Radiopharm';
        if injdosis <= 0
            error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
        end
    elseif isfield(info,'Private_0009_1038')
        injdosis = 1000000*info.Private_0009_1038;
        display='Private';
        if injdosis <= 0
            error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
        end
    else
        error('no RadionuclideTotalDose specified in dcm')
    end
    
    % write dcm info into tabel to check missing data
    table{1+i,1} = subject;
    table{1+i,1+1} = bodyweight;
    table{1+i,1+2} = injdosis;
    table{1+i,1+3} = display;
end

function [ptname,ptID,bodyweight, injdosis, acqtime, injtime, halftime, filelist_SUV]=LTNP_dcm2SUV(dicomdir,outputdir,outname)

    
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
    
%         if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose')
%             injdosis = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % in Bq
%             if injdosis <= 0 || injdosis > 400000000
%                 error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
%             end
%         elseif isfield(info,'Private_0009_1038')
%             injdosis = info.Private_0009_1038;
%             if injdosis <= 0 || injdosis > 400
%                 error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
%             end
%             injdosis = 1000000*info.Private_0009_1038;
%         elseif isfield(info,'Private_0009_103a')
%             injdosis = info.Private_0009_103a;
%             if injdosis <= 0 || injdosis > 400
%                 error('RadionuclideTotalDose specified in dcm is equal or lower than 0')
%             end
%             injdosis = 1000000*info.Private_0009_103a;
%         else
%             error('no RadionuclideTotalDose specified in dcm')
%         end
    
    
    
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
        
    dt=datetime(acqtime,'InputFormat','HHmmss')-datetime(injtime,'InputFormat','HHmmss.SSS');
    dt=seconds(dt); % in seconds
    
    if dt<3000
        error('difference between acquisition time and injection time is doubtfully small')
    end
    
    %RescaleIntercept=info.(dicomlookup('0028','1052'));
    %RescaleSlope= info.(dicomlookup('0028','1053'));
        
    % dcm2nii
    %command = strcat('/usr/local/bin/dcm2niix -o "',outputdir,'" -f "',outname,'" -p y -z n "',dicomdir,'"');
    % eval(['dcm2nii.exe -g N -m N -n Y -r Y -v N -x Y % -o',Outpudir,'',Inputfilename]) % Crop image
    %system(command)
    
    % Split niix into niiframes
    %tmp='TMP';
    cd(outputdir);
    mkdir('tmp')
    outputdir_TMP=fullfile(outputdir,'tmp');
    %nii_tot=char(fullfile(outputdir,[outname '.nii']));
    %spm_file_split(nii_tot,outputdir_TMP);
    %delete(nii_tot)
    %delete('*.json') 
    
    % Run dicom converter from spm
    nr_dcmframes=length(list);
    dcmframes = cell(nr_dcmframes,1);
    for i = 1:nr_dcmframes
       dcmframes{i,1} = fullfile(list(i).folder,list(i).name);
    end        

    % Initialise spm_jobman
    spm_jobman('initcfg')

    % Write job
    matlabbatch{1}.spm.util.import.dicom.data = dcmframes;
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = {outputdir_TMP};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.meta = 1;
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

    % Save batch
    cd(outputdir)
    save batch_dcm2nii matlabbatch
    
    % Run batchfile
    spm_jobman('run',matlabbatch)
       
    % calculate SUV image
    cd(outputdir_TMP)
    niiframes=dir(fullfile(outputdir_TMP,'*.nii'));
    nr_niiframes=length(niiframes);
    [~,Vref]=LCN12_read_image(niiframes(1).name);
    SUV = zeros(Vref.dim(1),Vref.dim(2),Vref.dim(3),nr_niiframes);
    for i = 1:nr_niiframes
        tmp = LCN12_read_image(niiframes(i).name,Vref);
        tmp = tmp * (bodyweight/injdosis) * (1/2)^(-dt/halftime);
        %tmp = (double(tmp)+RescaleIntercept)* RescaleSlope * (bodyweight/injdosis);
        SUV(:,:,:,i) = tmp;
        % Save Image
        %LCN12_write_image(tmp,fullfile(outputdir,[niiframes(i).name]),'SUV',Vref.dt(1),Vref); 
        LCN12_write_image(tmp,fullfile(outputdir,[outname '_0000' num2str(i) '.nii']),'SUV',Vref.dt(1),Vref);
        delete(niiframes(i).name)
    end
    delete('*.json')
    
    % Bring the filelist back
    filelist_SUV= dir(fullfile(outputdir,[outname '_0000*.nii']));
    
    % Delete temporary folder 
    cd(outputdir)
    rmdir('tmp','s')
    
end
function [ptname, ptID]=LTNP_dcm2niix(input_folder,output_folder,output_name)


% /this/is/my/input_folder
% /this/is/my/output_folder
% this_is_my_output_name_without_extension
    
% Loading dicom header
list=dir(input_folder);
list=list(~startsWith({list.name}, '.')); % removing hidden files from the list
dicomimg=fullfile(input_folder,list(10).name); % grabs tenth image (pseudorandom)
info = dicominfo(dicomimg);

% Get name, if inside dcm header
if isfield(info,'PatientName')
    if isfield(info.PatientName,'FamilyName') && isfield(info.PatientName,'GivenName')
        ptname=[info.PatientName.FamilyName '_' info.PatientName.GivenName];
    else
        ptname='unknown, please refer to ptID';
    end          
else
    ptname='unknown, please refer to ptID';
end

% Get ID, if inside dcm header
if isfield(info,'PatientID')
    ptID=info.PatientID;
else
    ptID='unknown, please refer to ptname';
end

% Run dcm2niix
command = strcat('/usr/local/bin/dcm2niix -o "',output_folder,'" -f "',output_name,'" -p y -z n "',input_folder,'"');
system(command)
% 
% % Run dicom converter from spm
% nr_frames=length(list);
% framelist = cell(nr_frames,1);
% for i = 1:nr_frames
%    framelist{i,1} = fullfile(list(i).folder,list(i).name);
% end        
% 
% % Initialise spm_jobman
% spm_jobman('initcfg')
% 
% % Write job
% matlabbatch{1}.spm.util.import.dicom.data = framelist;
% matlabbatch{1}.spm.util.import.dicom.root = 'flat';
% matlabbatch{1}.spm.util.import.dicom.outdir = {output_folder};
% matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
% matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
% matlabbatch{1}.spm.util.import.dicom.convopts.meta = 1;
% matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
% 
% % Save batch
% cd(output_folder)
% save batch_dcm2nii matlabbatch
% 
% % Run batchfile
% spm_jobman('run',matlabbatch)
% 
% % Rename
% niftis=dir(fullfile(output_folder,'*.nii'));
% 
% if length(niftis)==1
%     movefile(niftis(1).name,[output_name '.nii'])
%     delete('*.json')
% else
%     error('more or less than 1 nifti file found in the ouput_folder')
% end
    
end


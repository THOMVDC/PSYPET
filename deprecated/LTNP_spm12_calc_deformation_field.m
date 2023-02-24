function [def,invdef]=LTNP_spm12_calc_deformation_field(T1_path,outfolder)

%% Background
%
% Calculate the deformation and inverse deformation field of the warping of
% a T1 towards MNI 
%
% Input:
%       T1_path = string, absolute path to the T1 you want to calculate the deformation field to mni from
%       outfolder = string, absolute path to the folder you want the output to be stored
%
% Output:
%       def = string, absolute path to the deformation field
%       invdef = string, absolute path to the inverse deformation field
%
% Author:
%       Thomas Vande Casteele, KU Leuven
% 
% Dependency:
%       SPM12, Wellcome Trust Centre for Neuroimaging, University College London

%% Processing
% Defining input paths
spm_dir     = which('spm');
spm_dir     = spm_dir(1:end-6);
spm_dir_tpm = fullfile(spm_dir,'tpm');
%cat_dir     = fullfile(spm_dir,'toolbox','cat12');
addpath(spm_dir);
%addpath(cat_dir);

% Define output paths
[~,T1_name,T1_ext]=fileparts(T1_path);
out_image=fullfile(outfolder,[T1_name T1_ext]);
if strcmp(T1_path,out_image)==0
    copyfile(T1_path,outfolder) 
end
def=fullfile(outfolder,['y_' T1_name '.nii']);
invdef=fullfile(outfolder,['iy_' T1_name '.nii']);

% Initialise spm_jobman
spm_jobman('initcfg')

% Write job
matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr([out_image ',1']);
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii,1'));
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii,2'));
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii,3'));
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii,4'));
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii,5'));
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii,6'));
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % writes def and invdef
                                          
% Save batch
cd(outfolder)
save batch_coregister matlabbatch

% Run batchfile
spm_jobman('run',matlabbatch)
                                          
end
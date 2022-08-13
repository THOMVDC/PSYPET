function [out_image_warped]=LTNP_spm12_warp(image,deformation_field,outfolder,voxelsize)


% Input:
%       Absolute path to image to be warped
%       Absolute path to the deformation field that SPM12 will use
%       
% Output:
%       SPM warped image by the deformation field
%
% Author: 
%       Thomas Vande Casteele

% Defining paths
spm_dir     = which('spm');
spm_dir     = spm_dir(1:end-6);
cat_dir     = fullfile(spm_dir,'toolbox','cat12');
addpath(spm_dir);
addpath(cat_dir);

% Isotropic or anisotropic 
if size(voxelsize)==1 % iso
    voxelsize=[voxelsize voxelsize voxelsize];
end

% Grab image name, copy it
[~,input_image_name,input_image_ext]=fileparts(image);
out_image=fullfile(outfolder,[input_image_name input_image_ext]);
if strcmp(image,out_image)==0
    copyfile(image,outfolder) 
end

% Initialise spm_jobman
spm_jobman('initcfg')

% Write batch
matlabbatch = {};
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformation_field);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr([out_image ',1']);
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                           78   76  85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxelsize;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

% Save batchfile
cd(outfolder)
save batch_warping_PET matlabbatch

% run batchfile
spm_jobman('run',matlabbatch)

% Return path of warped image
out_image_warped=fullfile(outfolder,['w' input_image_name input_image_ext]);

end

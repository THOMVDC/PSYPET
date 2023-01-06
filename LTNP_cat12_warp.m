function [out_image_warped]=LTNP_cat12_warp(image,deformation_field,outfolder)

%% Background

% Warps an image with a deformation field by 4th degree B-spline interpolation 
%
% Input:
%       image = string, absolute path to image to be warped
%       deformation_field = string, absolute path to the deformation field
%
% Output:
%       out_image_warped = string, absolute path to the CAT12 warped image by the deformation field
%
% Author: 
%       Thomas Vande Casteele, KU Leuven
%       Dependency: CAT12, Friedrich  Schiller University  Jena,  Jena,  Germany

%% Processing

% Defining paths
spm_dir     = which('spm');
spm_dir     = spm_dir(1:end-6);
cat_dir     = fullfile(spm_dir,'toolbox','cat12');
addpath(spm_dir);
addpath(cat_dir);

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
matlabbatch{1}.spm.tools.cat.tools.defs.field1 = cellstr(deformation_field);
matlabbatch{1}.spm.tools.cat.tools.defs.images = cellstr([out_image ',1']);
matlabbatch{1}.spm.tools.cat.tools.defs.interp = 4;
matlabbatch{1}.spm.tools.cat.tools.defs.modulate = 0;

% Save batchfile
cd(outfolder)
save batch_warping_PET matlabbatch

% run batchfile
spm_jobman('run',matlabbatch)

% Return path of warped image
out_image_warped=fullfile(outfolder,['w' input_image_name input_image_ext]);

end
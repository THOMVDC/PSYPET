function [dst_image,rother_image]=LTNP_spm12_coregister(ref_image,mov_image,out_folder,other_image)

% Input:
%       Absolute path to anatomical T1 image
%       Absolute path to the output_folder of choice
%       
% Output:
%       SPM Registration of mov_image to ref_image, the prefix 'r' is added
%       to the filename
%
% Author: 
%       Thomas Vande Casteele

if nargin <4
    rother_image='';
elseif ischar(other_image)
    [~, other_image_name, other_image_ext]=fileparts(other_image);
    rother_image=fullfile(out_folder,['r' other_image_name other_image_ext]);
    copyfile(other_image,rother_image) 
else
    rother_image=cell(size(other_image));
    for i=1:length(other_image)
        [~, other_image_name, other_image_ext]=fileparts(other_image{i});
        rother_image{i}=fullfile(out_folder,['r' other_image_name other_image_ext]);
        copyfile(other_image{i},rother_image{i}) 
        rother_image{i}=[rother_image{i} ',1'];
    end
end

% Defining paths
spm_dir     = which('spm');
spm_dir     = spm_dir(1:end-6);
cat_dir     = fullfile(spm_dir,'toolbox','cat12');
addpath(spm_dir);
addpath(cat_dir);

% Copy the mov_image to the destination file
[~, mov_image_name, mov_image_ext]=fileparts(mov_image);
dst_image=fullfile(out_folder,['r' mov_image_name mov_image_ext]);
copyfile(mov_image,dst_image)

% Initialise spm_jobman
spm_jobman('initcfg')

% Write batchfile with ref_image, source_image and filelist
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(ref_image);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(dst_image); % is your moving image
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(rother_image);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Save batch
cd(out_folder)
save batch_coregister matlabbatch

% Run batchfile
spm_jobman('run',matlabbatch)



end
function [moved_image,moved_other_image]=LTNP_spm12_coregister_reslice_nn(ref_image,mov_image,out_folder,other_image)

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
    moved_other_image='';
elseif ischar(other_image)
    [~, other_image_name, other_image_ext]=fileparts(other_image);
    rother_image=fullfile(out_folder,['r' other_image_name other_image_ext]);
    copyfile(other_image,rother_image) 
    moved_other_image=fullfile(out_folder,['rr' other_image_name other_image_ext]); % define output name (r prefix added by spm after reslicing)
else
    rother_image=cell(size(other_image));
    moved_other_image=cell(size(other_image));
    for i=1:length(other_image)
        [~, other_image_name, other_image_ext]=fileparts(other_image{i});
        rother_image{i}=fullfile(out_folder,['r' other_image_name other_image_ext]);
        copyfile(other_image{i},rother_image{i}) 
        rother_image{i}=[rother_image{i} ',1'];
        moved_other_image{i}=fullfile(out_folder,['rr' other_image_name other_image_ext]);  % define output name (r prefix added by spm after reslicing)
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
moved_image=fullfile(out_folder,['rr' mov_image_name mov_image_ext]);  % define output name (r prefix added by spm after reslicing)

% Initialise spm_jobman
spm_jobman('initcfg')

% Write batchfile with ref_image, source_image and filelist
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(ref_image);
matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(dst_image); % is your moving image
matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(rother_image);
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0; % nearest neighbour
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Save batch
cd(out_folder)
save batch_coregister matlabbatch

% Run batchfile
spm_jobman('run',matlabbatch)

% % r prefix added to dst_image rother
% dst_image=fullfile(out_folder,['rr' mov_image_name mov_image_ext]); % RegisterReslice
% if ischar(other_image)
%     [~, other_image_name, other_image_ext]=fileparts(other_image);
%     rother_image=fullfile(out_folder,['rr' other_image_name other_image_ext]);
% elseif iscell(other_image)
%     for i=1:length(other_image)
%         [~, other_image_name, other_image_ext]=fileparts(other_image{i});
%         rother_image{i}=fullfile(out_folder,['rr' other_image_name other_image_ext]);
%     end
% end



end
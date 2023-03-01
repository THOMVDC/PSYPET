function [def,invdef]=LTNP_cat12_calc_deformation_field(T1_path,outfolder,voxelsize)

%% Background

% Calculates the deformation and inverse deformation field of the warping of
% a T1 towards MNI 
%
% Input:
%       T1_path = string, absolute path to the T1 you want to calculate the deformation field to mni from
%       outfolder = string, absolute path to the folder you want the output to be stored
% 
% Optional input:
%       voxelsize = numeric, isotropic voxelsize of the MNI space you want to warp to calculate the deformation field
%
% Output:
%       def = string, absolute path to the deformation field
%       invdef = string, absolute path to the inverse deformation field
%
% Author:
%       Thomas Vande Casteele, KU Leuven
% 
% Dependency:
%       CAT12, Friedrich  Schiller University  Jena,  Jena,  Germany

%% Processing
% Defining input paths
spm_dir     = which('spm');
spm_dir     = spm_dir(1:end-6);
spm_dir_tpm = fullfile(spm_dir,'tpm');
cat_dir     = fullfile(spm_dir,'toolbox','cat12');
addpath(spm_dir);
addpath(cat_dir);

% Define output paths
[~,T1_name,T1_ext]=fileparts(T1_path);
out_image=fullfile(outfolder,[T1_name T1_ext]);
if strcmp(T1_path,out_image)==0
    copyfile(T1_path,outfolder) 
end
def=fullfile(outfolder,'mri',['y_' T1_name '.nii']);
invdef=fullfile(outfolder,'mri',['iy_' T1_name '.nii']);

% Voxelsize
if nargin<3
    voxelsize=1.5;
end

% Initialise spm_jobman
spm_jobman('initcfg')

% Write job
matlabbatch{1}.spm.tools.cat.estwrite.data = cellstr([out_image ',1']);
%matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0; % if you want parallelization, set this to the number of cpu's you want to dedicate
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
%matlabbatch{1}.spm.tools.cat.estwrite.opts.ngaus = [1 1 2 3 4 2];
%matlabbatch{1}.spm.tools.cat.estwrite.opts.warpreg = [0 0.001 0.5 0.05 0.2];
%matlabbatch{1}.spm.tools.cat.estwrite.opts.bias.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.opts.acc.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.opts.redspmres = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.fixed = [1 0.1];
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.T1 = cellstr(fullfile(spm_dir,'toolbox','FieldMap','T1.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask =  cellstr(fullfile(spm_dir,'toolbox','FieldMap','brainmask.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas = cellstr(fullfile(cat_dir,'templates_1.50mm','cat.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = cellstr(fullfile(cat_dir,'templates_1.50mm','Template_1_IXI555_MNI152.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm =cellstr(fullfile(cat_dir,'templates_1.50mm','Template_1_IXI555_MNI152.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = cellstr(fullfile(cat_dir,'templates_1.50mm','Template_0_IXI555_MNI152_GS.nii'));   
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = voxelsize;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
                                          
% Save batch
cd(outfolder)
save batch_coregister matlabbatch

% Run batchfile
spm_jobman('run',matlabbatch)
                                          
end

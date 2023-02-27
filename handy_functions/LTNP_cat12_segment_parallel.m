function LTNP_cat12_segment_parallel(T1, output_folder, voxelsize, wmhc, atlas)

%% Background
% 
% Segmentats multiple T1s with CAT12, results can be found in the same folder
% as the respective T1 images
%
% Input:
%       T1 = cell(n,1) with the n absolute paths to anatomical T1 images in nifti
%       output_folder = string, absolute path to the output_folder of choice
%       wmhc = numeric, should be equal to 1 (wmh map output) or to 0 (no wmh map output)
%              ! the cat.extopts.WMHC should be set to 2 (considers WMH and corrects segmentation to WM)
%       voxelsize = numeric, voxelsize of the space you want the normalised output to be in
%
% Optional input:
%       atlas= string, should be [] | 'cobra' | 'hammers' | 'neuromorphometrics' | 'lbpa40' |
%       '/path/to/own_atlas'
%
% Author: 
%       Thomas Vande Casteele, KU Leuven
%
% Dependency:
%       CAT12, Friedrich  Schiller University  Jena,  Jena,  Germany

%% Processing
% Defining paths
spm_dir     = which('spm');
spm_dir     = spm_dir(1:end-6);
spm_dir_tpm = fullfile(spm_dir,'tpm');
cat_dir     = fullfile(spm_dir,'toolbox','cat12');
addpath(spm_dir);
addpath(cat_dir);
cat12('expert');

% Making output folder with a copy of T1
if ischar(T1)
    [~,image_name,image_ext] = fileparts(T1);
    T1_copy=fullfile(output_folder,[image_name image_ext]);
    if strcmp(T1,T1_copy)==0
        copyfile(T1,output_folder)
    end
    T1_copy={[T1_copy ',1']};
else
    T1_copy=cell(size(T1));
    for i=1:length(T1)
        [~, image_name, image_ext]=fileparts(T1{i});
        T1_copy{i}=fullfile(output_folder,[image_name image_ext]);
        if strcmp(T1{i},T1_copy{i})==0
            copyfile(T1{i},output_folder)
        end
        T1_copy{i}=[T1_copy{i} ',1'];
    end
end

% Check if and which atlas is specified
if nargin > 4 % then atlas is specified
    neuromorphometrics=0;
    lpba40=0;
    cobra=0;
    hammers=0;
    isbr=0;
    aal=0;
    mori=0;
    anatomy=0;
    if strcmp(atlas,'neuromorphometrics')
        neuromorphometrics=1;
    elseif strcmp(atlas,'lpba40')
        lpba40=1;
    elseif strcmp(atlas,'cobra')
        cobra=1;
    elseif strcmp(atlas,'isbr')
        isbr=1;
    elseif strcmp(atlas,'aal3')
        aal=1;
    elseif strcmp(atlas,'mori')
        mori=1;
    elseif strcmp(atlas,'anatomy3')
        anatomy=1;
    end
else % no atlas specified, set a default mode
    neuromorphometrics=1;
    lpba40=0;
    cobra=1;
    hammers=1;
    isbr=0;
    aal=0;
    mori=0;
    anatomy=0;
end

% Set parameters for white matter hyperintensities map output
if wmhc==0
    WMHC_seg=1;
    WMH_native = 0;
    WMH_warped = 0;
    WMH_mod = 0;
	WMH_dartel = 0;
elseif wmhc==1
    WMHC_seg=3;
	WMH_native = 1;
    WMH_warped = 1;
    WMH_mod = 1;
	WMH_dartel = 1;
end

% Initialise spm_jobman
spm_jobman('initcfg')

% Write batchfile CAT12.6
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.expertgui = 2; % force
%expert mode
matlabbatch{1}.spm.tools.cat.estwrite.data = T1_copy;
%matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 8; % if you want parallelization, set this to the number of cpu's you want to dedicate
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
%matlabbatch{1}.spm.tools.cat.estwrite.opts.ngaus = [1 1 2 3 4 2];
%matlabbatch{1}.spm.tools.cat.estwrite.opts.warpreg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.tools.cat.estwrite.opts.bias.biasstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.acc.accstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.opts.redspmres = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = WMHC_seg;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.fixed = [1 0.1];
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.T1 = cellstr(fullfile(spm_dir,'toolbox','FieldMap','T1.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask =  cellstr(fullfile(spm_dir,'toolbox','FieldMap','brainmask.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas = cellstr(fullfile(cat_dir,'templates_1.50mm','cat.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = cellstr(fullfile(cat_dir,'templates_1.50mm','Template_1_IXI555_MNI152.nii'));
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm =cellstr(fullfile(cat_dir,'templates_1.50mm','Template_1_IXI555_MNI152.nii'));
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
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = neuromorphometrics;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = lpba40;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = cobra;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = hammers;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = isbr;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal = aal;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = mori;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy = anatomy;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = WMH_native;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = WMH_warped;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = WMH_mod;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = WMH_dartel;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];


% Save batchfile
cd(output_folder)
save batch_segment matlabbatch

% Run batchfile
spm_jobman('run',matlabbatch)

% Load CAT12 measures
%CAT12_vol_thick=load(fullfile(output_folder,'report',['cat_' image_name '.mat']));

% Catch TIV
%str=CAT12_vol_thick.S.catlog{end-4, 1};
%IQR=extractBetween(str,'(IQR):  ','%');
%thickness=cell2mat(CAT12_vol_thick.S.subjectmeasures.dist_thickness);
%TIV=CAT12_vol_thick.S.subjectmeasures.vol_TIV;

% Grab path of deformation fields
%deformation_field = fullfile(output_folder,'mri',['y_' image_name image_ext]);
%invdeformation_field = fullfile(output_folder,'mri',['iy_' image_name image_ext]);

% TIV=cell(size(T1));
% deformation_field=cell(size(T1));
% invdeformation_field=cell(size(T1));
% for t=1:length(T1)
%     [~, image_name, image_ext]=fileparts(T1{i});
%     CAT12_vol_thick=load(fullfile(output_folder,'report',['cat_' image_name '.mat']));
%     TIV{t}=CAT12_vol_thick.S.subjectmeasures.vol_TIV;
%     deformation_field{t} = fullfile(output_folder,'mri',['y_' image_name image_ext]);
%     invdeformation_field{t} = fullfile(output_folder,'mri',['iy_' image_name image_ext]);
% end


end

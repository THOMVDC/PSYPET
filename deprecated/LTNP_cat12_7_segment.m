function [IQR, TIV, thickness, deformation_field, invdeformation_field]=LTNP_cat12_7_segment(T1, output_folder, voxelsize, wmhc, atlas)

%% Background
% 
% Segmentents a specified T1 with CAT12, results can be found in the same folder
% as the T1 image
%
% Input:
%       T1 = string, absolute path to anatomical T1 images in nifti
%       output_folder= string, absolute path, to the output_folder of choice
%       wmhc = numeric, should be equal to 1 (wmh map output) or to 0 (no wmh map output)
%              ! the cat.extopts.WMHC should be set to 2 (considers WMH and corrects segmentation to WM)
%       voxelsize = numeric, voxelsize of the space you want the normalised output to be in
%
% Optional input:
%       atlas= string, should be [] | 'cobra' | 'hammers' | 'neuromorphometrics' | 'lbpa40' |
%       '/path/to/own_atlas'
% Output:
%       IQR = numeric, overall quality rating of T1 provided by CAT12
%       TIV = numeric, estimated TIV provided by CAT12
%       thickness = numeric, estimated average cortical thickness as
%                   provided by CAT12
%       deformation_field = string, absolute path to the deformation field
%       invdeformation_field = string, absolute path to the inverse deformation field
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

% Making output folder with a copy of T1
[~,image_name,image_ext] = fileparts(T1);
T1_copy=fullfile(output_folder,[image_name image_ext]);
if strcmp(T1,T1_copy)==0
    copyfile(T1,output_folder)
end

% Check if and which atlas is specified
if nargin > 4 % then atlas is specified
    own_atlas={''};
    neuromorphometrics=0;
    lpba40=0;
    cobra=0;
    hammers=0;
    isbr=0;
    aal3=0;
    mori=0;
    anatomy3=0;
    julichbrain=0;
    schaefer100=0;
    schaefer200=0;
    schaefer400=0;
    schaefer600=0;
    if strcmp(atlas,'neuromorphometrics')
        neuromorphometrics=1;
    elseif strcmp(atlas,'lpba40')
        lpba40=1;
    elseif strcmp(atlas,'cobra')
        cobra=1;
    elseif strcmp(atlas,'isbr')
        isbr=1;
    elseif strcmp(atlas,'aal3')
        aal3=1;
    elseif strcmp(atlas,'mori')
        mori=1;
    elseif strcmp(atlas,'anatomy3')
        anatomy3=1;
    elseif strcmp(atlas,'julichbrain')
        julichbrain=1;
    elseif strcmp(atlas,'schaefer100')
        schaefer100=1;
    elseif strcmp(atlas,'schaefer200')
        schaefer200=1;
    elseif strcmp(atlas,'schaefer400')
        schaefer400=1;
    elseif strcmp(atlas,'schaefer600')
        schaefer600=1;
    else
        own_atlas=cellstr([atlas ',1']);
    end
else % no atlas specified, set a default mode
    own_atlas={''};
    neuromorphometrics=1;
    lpba40=0;
    cobra=1;
    hammers=1;
    isbr=0;
    aal3=0;
    mori=0;
    anatomy3=0;
    julichbrain=0;
    schaefer100=0;
    schaefer200=0;
    schaefer400=0;
    schaefer600=0;
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

% Write batchfile CAT12.7
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.expertgui = 1; % force
%expert mode
matlabbatch{1}.spm.tools.cat.estwrite.data =  cellstr([T1_copy ',1']);
matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 8;
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = cellstr(fullfile(spm_dir_tpm,'TPM.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = WMHC_seg;% set white matter hypertensities as a separate class (WMHC=0 no correction, WMHC=1 set as WM temporarly for processing, WMHC=2 set as WM permanently, WMHC=3 set as own class
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = cellstr(fullfile(cat_dir,'templates_volumes','Template_0_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = voxelsize;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbt2x';
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtlas = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.collcorr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist = 1.33333333333333;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = neuromorphometrics;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = lpba40;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = cobra;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = hammers;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = isbr;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = aal3;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = mori;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3 = anatomy3;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = julichbrain;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = schaefer100;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = schaefer200;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = schaefer400;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = schaefer600;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = own_atlas;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = WMH_native;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = WMH_warped;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = WMH_mod;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = WMH_dartel;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 1;

% Save batchfile
cd(output_folder)
save batch_segment matlabbatch

% Run batchfile
spm_jobman('run',matlabbatch)

% Load CAT12 measures
CAT12_vol_thick=load(fullfile(output_folder,'report',['cat_' image_name '.mat']));

% Catch IQR (overall quality rating), thickness and TIV
str=CAT12_vol_thick.S.catlog{end-4, 1};
IQR=extractBetween(str,'(IQR):  ','%');
thickness=cell2mat(CAT12_vol_thick.S.subjectmeasures.dist_thickness);
TIV=CAT12_vol_thick.S.subjectmeasures.vol_TIV;

% Grab path of deformation fields
deformation_field = fullfile(output_folder,'mri',['y_' image_name image_ext]);
invdeformation_field = fullfile(output_folder,'mri',['iy_' image_name image_ext]);

end
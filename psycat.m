function [segmentation,invdef]=psycat(T1,atlas,outfolder)

% Set defaults
WMHC=1;
voxelsize=1.5; %isotropic

% Get file name of T1
[~, T1name, T1ext]=fileparts(T1);

% Segment T1 and grab the inverse deformation field (mni -> patient space)
[~,~,invdef]=LTNP_cat12_6_segment(T1, outfolder, voxelsize, WMHC, atlas); % also creates the deformation field (warping parameters) to MNI

% Grab tissue probability maps
GM_path=fullfile(outfolder,'mri',['p1' T1name T1ext]);
WM_path=fullfile(outfolder,'mri',['p2' T1name T1ext]);
CSF_path=fullfile(outfolder,'mri',['p3' T1name T1ext]);

% Make tissue masks
[GMmask_path, WMmask_path, CSFmask_path,~,~]=LTNP_make_labelimage_v4(GM_path,WM_path,CSF_path,outfolder);

% Make final atlas based segmentation
atlas_path=fullfile(outfolder,'mri_atlas',[atlas '_' T1name T1ext]);
[~,~,~,segmentation]=LTNP_create_rbv_segm_v4(GMmask_path, WMmask_path, CSFmask_path, atlas_path, outfolder);
   
end

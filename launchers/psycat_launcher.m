%% Background
% Preprocess T1 dcm folder independently from psypet

%% Settings
dcmdir='';
outfolder='';
WMHC=1;
voxelsize=1.5; %isotropic
subj='B002';

%% Processing
% Give a name to the (to be) processed image
T1name=['T1_' subj '.nii'];
outfolder_t1=fullfile(outfolder,subj);
cd(outfolder);
mkdir(subj);

% Convert to nifti, center
[accT1nii, ~, ~]=LTNP_preproc_T1_spm(dcmdir,outfolder_t1,T1name); 

% delete qform
LTNP_delete_qform(accT1nii,accT1nii) 

% Segment T1 and grab the inverse deformation field (mni -> patient space)
[~,~,~,GM_path, WM_path, CSF_path]=LTNP_cat12_segment(T1, outfolder_t1, voxelsize, WMHC, atlas); % also creates the deformation field (warping parameters) to MNI

% Make tissue masks
[GMmask_path, WMmask_path, CSFmask_path,~,~]=LTNP_make_labelimage_v4(GM_path,WM_path,CSF_path,outfolder_t1);

% Make final atlas based segmentation
atlas_path=fullfile(outfolder,'mri_atlas',[atlas '_' T1name T1ext]);
[~,~,~,segmentation]=LTNP_create_rbv_segm_v4(GMmask_path, WMmask_path, CSFmask_path, atlas_path, outfolder_t1);
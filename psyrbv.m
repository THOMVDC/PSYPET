function [rrrSUVnii_pvc_RBV_path,wrrrSUVnii_pvc_RBV_path]=psyrbv(T1,subjectcode,subjdir_ANAT,rrrSUVnii,subjdir_PET)

% Define rbv_script_dir (should be located in the same path as here)
[rbv_script_dir,~,~]=fileparts(mfilename('fullpath'));

% Define directories
subjdir_ANATsegm = fullfile(subjdir_ANAT,'CAT12');
subjdir_ANATmask = fullfile(subjdir_ANAT,'MASKS2');
subjdir_ANATatlas= fullfile(subjdir_ANAT,'RBV_ATLASSES2');
subjdir_PETpvcrbv = fullfile(subjdir_PET,'PVC_RBV2');
mkdir(subjdir_PETpvcrbv);
mkdir(subjdir_ANATmask);
mkdir(subjdir_ANATatlas);
subjdir_PETw = fullfile(subjdir_PET,'WARPED');
 
% Grab name of T1
[~,ref_image_name,ref_image_ext]=fileparts(T1);

% Segment
GM_path=fullfile(subjdir_ANATsegm,'mri',['p1' ref_image_name ref_image_ext]); % GM in patientspace
WM_path=fullfile(subjdir_ANATsegm,'mri',['p2' ref_image_name ref_image_ext]); % WM in patientspace
CSF_path=fullfile(subjdir_ANATsegm,'mri',['p3' ref_image_name ref_image_ext]); % CSF in patientspace
deformation_field=fullfile(subjdir_ANATsegm,'mri',['y_' ref_image_name ref_image_ext]); % CSF in patientspace
%WMH_path=fullfile(subjdir_ANATsegm,'mri',['p7' ref_image_name ref_image_ext]); % W
wGM_path=fullfile(subjdir_ANATsegm,'mri',['mwp1' ref_image_name ref_image_ext]); % GM in patientspace
wWM_path=fullfile(subjdir_ANATsegm,'mri',['mwp2' ref_image_name ref_image_ext]); % WM in patientspace
wCSF_path=fullfile(subjdir_ANATsegm,'mri',['mwp3' ref_image_name ref_image_ext]); % CSF in patientspace
%wWMH_path=fullfile(subjdir_ANATsegm,'mri',['mwp7' ref_image_name ref_image_ext]); % W

% Make new masks
[GMmask_path, WMmask_path, CSFmask_path,~,~]=LTNP_make_labelimage(GM_path,WM_path,CSF_path,subjdir_ANATmask);
[~,~,~,~,~]=LTNP_make_labelimage(wGM_path,wWM_path,wCSF_path,subjdir_ANATmask);

% Grab output atlas generated from segmentation (you'll need it later)
atlas_mni_name='neuromorphometrics';
path_to_atlasses_cat12=fullfile(subjdir_ANATsegm,'mri_atlas'); % Grab atlasses in patient space
atlas_path=fullfile(path_to_atlasses_cat12,[atlas_mni_name '_' ref_image_name ref_image_ext]);
%watlas_path=fullfile(path_to_atlasses_cat12,['w' atlas_mni_name '_' ref_image_name ref_image_ext]);

% Create rbv segmentation and brain mask
[rbv_atlas_path,~,~]=LTNP_create_rbv_segm_v13(GMmask_path, WMmask_path, CSFmask_path, atlas_path, subjdir_ANATatlas);
copyfile(rbv_atlas_path,subjdir_PETpvcrbv)

% Bring RBV atlas to mni
LTNP_spm12_warp(rbv_atlas_path,deformation_field,subjdir_PETpvcrbv);

% Execute RBV
cd(subjdir_PETpvcrbv) % change directory to output directory for RBV script
[rrrSUVnii_pvc_RBV_path]=LTNP_PVC_RBV(subjectcode,rbv_script_dir,rrrSUVnii,rbv_atlas_path,subjdir_PETpvcrbv);

% Move files to output_folder
mov=dir([subjdir_PET '/PVC*.*']);
for d=1:numel(mov)
    movefile(fullfile(subjdir_PET,mov(d).name),subjdir_PETpvcrbv);
end
[~,rrrSUVnii_pvc_RBV_path_img,rrrSUVnii_pvc_RBV_path_ext]=fileparts(rrrSUVnii_pvc_RBV_path);
rrrSUVnii_pvc_RBV_path=fullfile(subjdir_PETpvcrbv,[rrrSUVnii_pvc_RBV_path_img rrrSUVnii_pvc_RBV_path_ext]);
wrrrSUVnii_pvc_RBV_path=LTNP_spm12_warp(rrrSUVnii_pvc_RBV_path,deformation_field,subjdir_PETw);

end
%function psypetref(T1_fs,T1_cat12, FLAIR, WM, WML, ref_VOI_mni)
% stats

function SUVR_path=make_SUVR_with_fs(accT1nii, rSUVnii_path, parcellation, ref_VOI_fs, output_folder)

%% Background
%
% accT1nii=brain.nii file from fastsurfer output
% rSUVnii_path=realigned and summated PET acquisition
% parcellation=aseg+aparc file from fastsurfer output
% ref_VOI=reference volume in patient space

%% 1/ Coregister SUV to T1
%       ----------------------

% Check if participant PET is corresponding to participant T1
% ptname_T1(isspace(ptname_T1)) = []; % remove spaces
% ptname_T1=lower(ptname_T1); % lower all letters
% ptname_T1=replace(ptname_T1,'_',''); % remove underscores
% ptname_PET(isspace(ptname_PET)) = []; % remove spaces
% ptname_PET=lower(ptname_PET); % lower all letters
% ptname_PET=replace(ptname_PET,'_',''); % remove underscores
% 
% if strcmp(ptname_T1,ptname_PET) || strcmp(ptID_T1,ptID_PET) || strcmp(ptname_T1,ptID_PET) || strcmp(ptname_PET,ptID_T1)
%    fprintf(fid,'T1 and PET corresponding');
% else
%    fprintf(fid,'T1 and PET not corresponding');
% end

% Talk to config file
%fprintf(fid,'\n');
%fprintf(fid,'Coregistration started \n');

% First center PET image (T1 already aligned on the AC)
LTNP_center(rSUVnii_path)

% Let the magic happen
%rigid='r';
rrrSUVnii=LTNP_spm12_coregister_reslice(accT1nii, rSUVnii_path, output_folder); % first rigid coreg with spm
%rrrSUVnii=LTNP_ANTs_coregister(accT1nii,rrSUVnii,rigid, subjdir_PETcor);  % now second coreg with ants

% Talk to config file
%fprintf(fid,'Coregistration ended \n');
%% 2/ Apply PVC SUV


% Talk to log
fprintf(fid,'Applying Region Based partial volume correction');
        
% Execute RBV
cd(subjdir_PETpvcrbv) % change directory to output directory for RBV script
[rrrSUVnii_pvc_RBV_path]=LTNP_PVC_RBV(subjectcode,rbv_script_dir,rrrSUVnii,parcellation,subjdir_PETpvcrbv);

% Move files to output_folder
mov=dir([subjdir_PET '/PVC*.*']);
for d=1:numel(mov)
    movefile(fullfile(subjdir_PET,mov(d).name),subjdir_PETpvcrbv);
end
[~,rrrSUVnii_pvc_RBV_path_img,rrrSUVnii_pvc_RBV_path_ext]=fileparts(rrrSUVnii_pvc_RBV_path);
rrrSUVnii_pvc_RBV_path=fullfile(subjdir_PETpvcrbv,[rrrSUVnii_pvc_RBV_path_img rrrSUVnii_pvc_RBV_path_ext]);

% Talk to log
fprintf(fid,'Region Based partial volume correction done');
%% 3/ Make SUVR image
 
% Read data
[SUV, Vref]=LCN12_read_image(rrrSUVnii_pvc_RBV_path);
%GMimg=LCN12_read_image(GMpath,Vref);
refVOIimg=LCN12_read_image(ref_VOI_fs,Vref);

% Threshold
ref_mask=refVOIimg>0.5;

% Use the LTNP_calculate_SUVR functions   
ref_value        = nanmean(SUV(ref_mask>0));

% Calculate SUVR
SUVR             = SUV./ref_value;

% Save SUVR, ref mask, rWMHimg_thresholded
[~, SUV_name, ~]=fileparts(rrrSUVnii_pvc_RBV_path);
SUVR_path = fullfile(out_folder,['SUVR_' SUV_name '_' ref_VOI_name '_' subjectcode '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min.nii']); 
LCN12_write_image(SUVR,SUVR_path,'SUVR',Vref.dt(1),Vref); 
     
end


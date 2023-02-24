% subjdir_ANATsegm='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/T1/B014/CAT12';
% ref_image_name='accT1_B014';
% ref_image_ext='.nii';
% 

WMHC=1;
 
% out_folder='/Users/tvdcas1/Desktop/mask4';
% 
% LTNP_make_labelimage(GM_path,WM_path,CSF_path,out_folder)
% 

ref_VOI_mni ='/Volumes/LaCie/Thomas/Projects/SCRIPTS/PSYPET/templates/Mask_VOI_SO_AtlasspaceSPM.nii'; %  TRACER DEPENDENT
[~,ref_VOI_name,~] = fileparts(ref_VOI_mni);
    
% tracer= 'UCBJ'; %  TRACER DEPENDENT
% atlas_mni='neuromorphometrics';
% atlas_csv='neuromorphometrics';
% 
% output_folder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/UCBJ_4s'; %  TRACER DEPENDENT
% %maindir_T1='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/T1';
maindir_T1='/Volumes/LaCie/Margot/aline/T1_cat12';
maindir_PET='/Volumes/LaCie/Margot/aline/Delva_UCBJ'; %  TRACER DEPENDENT
subjects=dir([maindir_T1 '/HV_0*']);
tracer= 'UCBJ'; %  TRACER DEPENDENT
GM_threshold                 = 0.3;
WM_threshold                 = 0.2;
SUVR_start_time='60';
SUVR_end_time='90';


% s=[35,40,58]

for s=1:length(subjects)
    
    % Grab subj
    subj=subjects(s).name;
    disp(subj)
    
    % Grab SUV_dir
    SUV_dir=dir(fullfile(maindir_PET,subj,'psypet*',tracer,'COREG',['rrrSUV_' tracer '_' subj '.nii'])); % change path construction to grab the desired SUVR
    
    % Check if there is only one PET file
    if length(SUV_dir)==1
        PET=fullfile(SUV_dir(1).folder,SUV_dir(1).name);
        [subjdir_PET,~,~]=fileparts(PET);
        [subjdir_PET,~,~]=fileparts(subjdir_PET);
        subjdir_PETsuvr=fullfile(subjdir_PET,'SUVR');
    else
        fprintf(subj,'0 or more than one SUVR image available \n');
    end
    subjdir_ANAT=fullfile(maindir_T1,subj);
    T1=fullfile(subjdir_ANAT,['accT1_' subj '.nii']);
    [~,ref_image_name,ref_image_ext]=fileparts(T1);
    subjdir_ANATsegm=fullfile(subjdir_ANAT,'CAT12');
    
    % Run rbv
    [rrrSUVnii_pvc_RBV_path,wrrrSUVnii_pvc_RBV_path]=psyrbv(T1,subj,subjdir_ANAT,PET,subjdir_PET);

    % Grab probability maps
    GM_path=fullfile(subjdir_ANATsegm,'mri',['p1' ref_image_name ref_image_ext]); % GM in patientspace
    WM_path=fullfile(subjdir_ANATsegm,'mri',['p2' ref_image_name ref_image_ext]); % WM in patientspace
    CSF_path=fullfile(subjdir_ANATsegm,'mri',['p3' ref_image_name ref_image_ext]); 
    WMH_path=fullfile(subjdir_ANATsegm,'mri',['p7' ref_image_name ref_image_ext]);
    wGM_path=fullfile(subjdir_ANATsegm,'mri',['p1' ref_image_name ref_image_ext]); % GM in mnispace
    wWM_path=fullfile(subjdir_ANATsegm,'mri',['p2' ref_image_name ref_image_ext]); % WM in mnispace
    wCSF_path=fullfile(subjdir_ANATsegm,'mri',['p3' ref_image_name ref_image_ext]); 
    wWMH_path=fullfile(subjdir_ANATsegm,'mri',['p7' ref_image_name ref_image_ext]);
    
    fprintf('%s : processing done, no error occured \n',subj);
    ref_VOI_patient=dir(fullfile(subjdir_PET,'SUVR',['ref_mask_SUVR_rrrSUV_UCBJ_' subj '_Mask_VOI_SO_AtlasspaceSPM_UCBJ_' subj '_*min_*min_refVOI.nii']));
    ref_VOI_mni=dir(fullfile(subjdir_PET,'SUVR',['ref_mask_SUVR_rrrSUV_UCBJ_' subj '_Mask_VOI_SO_AtlasspaceSPM_UCBJ_' subj '_*min_*min_refVOI.nii']));
    ref_VOI_patient=fullfile(ref_VOI_patient(1).folder,ref_VOI_patient(1).name);
    ref_VOI_mni=fullfile(ref_VOI_mni(1).folder,ref_VOI_mni(1).name);
    BRAINmask_patient=fullfile(maindir_T1,subj,'MASKS2',['GMWMCSF_mask_paccT1_' subj '.nii']);
    BRAINmask_mni=fullfile(maindir_T1,subj,'MASKS2',['GMWMCSF_mask_mwpaccT1_' subj '.nii']);
    
    %Make SUVR
    PETMR_outdata=cell(2,3);
    PETMR_outdata(1,1)={rrrSUVnii_pvc_RBV_path};
    PETMR_outdata(2,1)={wrrrSUVnii_pvc_RBV_path};
    for k=1:2
        [SUV, Vref]=LCN12_read_image(PETMR_outdata{k,1});
        if k==1 % we'll read patient space data
            GMimg=LCN12_read_image(GM_path,Vref);
            WMimg=LCN12_read_image(WM_path,Vref);
            WMH=WMH_path;
            refVOIimg=LCN12_read_image(ref_VOI_patient,Vref);
            brain_mask=LCN12_read_image(BRAINmask_patient,Vref); 
        else % we'll read mni space data
            GMimg=LCN12_read_image(wGM_path,Vref);
            WMimg=LCN12_read_image(wWM_path,Vref);
            WMH=wWMH_path;
            refVOIimg=LCN12_read_image(ref_VOI_mni,Vref);
            brain_mask=LCN12_read_image(BRAINmask_mni,Vref);
            %atlas_path=watlas_path;
        end
        if WMHC==1
            WMHimg=LCN12_read_image(WMH,Vref);
        else
            WMHimg=[];
        end

        % Threshold
        refVOI_thresholded=refVOIimg>0.5;
        brain_mask_thresholded=brain_mask>0.5;

        % Use the LTNP_calculate_SUVR functions   
        if strcmp(tracer,'FLUT')
            [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_flut(SUV,refVOI_thresholded,GMimg,GM_threshold,brain_mask_thresholded);  
        elseif strcmp(tracer,'UCBJ')
            [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_ucbj(SUV,refVOI_thresholded,WMimg,WM_threshold,brain_mask_thresholded,WMHimg);
        elseif strcmp(tracer,'MK62')
            [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_mk62(SUV,refVOI_thresholded,GMimg,GM_threshold,brain_mask_thresholded);
        elseif strcmp(tracer,'PIB')
            [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_pib(SUV,refVOI_thresholded,GMimg,GM_threshold,brain_mask_thresholded);
        else
            [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR(PET,refVOI_thresholded,brain_mask_thresholded);
        end

        % Save SUVR, ref mask, nr_voxels_refVOI and ref_value
        PETMR_outdata{k,2}=nr_voxels_refVOI;
        PETMR_outdata{k,3}=ref_value;
        [~, SUV_name, ~]=fileparts(PETMR_outdata{k,1});
        outputname_refmask = fullfile(subjdir_PETsuvr,['ref_mask_SUVR_' SUV_name '_' ref_VOI_name '_' tracer '_' subj '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min_refVOI.nii']);
        SUVR_path = fullfile(subjdir_PETsuvr,['SUVR_' SUV_name '_' ref_VOI_name '_' tracer '_' subj '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min.nii']); 
        LCN12_write_image(SUVR_img,SUVR_path,'SUVR',Vref.dt(1),Vref); 
        LCN12_write_image(ref_mask,outputname_refmask,'SUVR_refVOI',Vref.dt(1),Vref);
    end
  
end

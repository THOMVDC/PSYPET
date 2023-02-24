WMHdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/WML/';
T1PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET2.0/UCBJ/';
FLAIRdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/FLAIR/';
BRAINmask_path='';
PSYPETFSdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPETFASTSURFER/';
FSdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/FASTSURFER/no_biais_correction/fastsurfer_output_1/';
rbv_script_dir='/Volumes/LaCie/Thomas/Projects/SCRIPTS/PSYPET';

subjects=dir([T1PETdir 'B0*']);
for s=1:length(subjects)
    
    % Grab subject code
    subj=subjects(s).name;
    
    % Define and make output folder
    in_folder=dir(fullfile(T1PETdir,subj,'psypet_v2.0_processed_CAT12.7_*2022'));
    mkdir(fullfile(in_folder.folder,in_folder.name),'FS');
    out_folder_FS=fullfile(in_folder.folder,in_folder.name,'FS');
    out_folder_suvr=fullfile(in_folder.folder,in_folder.name,'UCBJ','SUVR');
    ref_VOI_name='brainstem';
        
    % Grab corresponding images
    rSUVnii_path=fullfile(in_folder.folder,in_folder.name,'UCBJ','REALIGNED',['rSUV_UCBJ_' subj '.nii']);
    %FLAIRpath=fullfile(FLAIRdir,['FLAIR_' subj '.nii']);
%     T1image_name=['accT1_' subj '.nii'];
%     T1path_cat12=fullfile(in_folder.folder,in_folder.name,'ANAT',T1image_name);
    T1path=fullfile(FSdir,subj,'mri','brain.nii');
    %PETpath=fullfile(in_folder.folder,in_folder.name,'UCBJ','COREG',['rrrSUV_UCBJ_' subj '.nii']);
    PETpath=fullfile(PSYPETFSdir,['PVC_coreg_SUV_aseg_' subj '_PET_PVC_RBV_65mm_in_seg.nii']);
%     GMpath=fullfile(in_folder.folder,in_folder.name,'ANAT','CAT12','mri',['p1' T1image_name]); % GM in patientspace
%     WMpath=fullfile(in_folder.folder,in_folder.name,'ANAT','CAT12','mri',['p2' T1image_name]); % WM in patientspace
%     CSFpath=fullfile(in_folder.folder,in_folder.name,'ANAT','CAT12','mri',['p3' T1image_name]); % CSF in patientspace
%     WMHpath=fullfile(WMHdir,['fs' subj '_lesions.nii']);
    %autoWMHpath=fullfile(in_folder.folder,in_folder.name,'ANAT','MASKS',['mask_p7accT1_' subj '.nii']);
    %WMHpath=fullfile(subjdir_ANATsegm,'mri',['p7' ref_image_name ref_image_ext]); % White matt
    %ref_VOI_patient=fullfile(in_folder.folder,in_folder.name,'UCBJ','WARPED','wMask_VOI_SO_AtlasspaceSPM.nii');
    ref_VOI_patient=fullfile(FSdir,subj,'mri',[ref_VOI_name '.nii']);
    %BRAINmask_path=fullfile(in_folder.folder,in_folder.name,'ANAT','MASKS',['GMWMCSF_mask_p' T1image_name]);
    %ref_mask_SO_without_autoWMH_path=fullfile(in_folder.folder,in_folder.name,'UCBJ','SUVR',['ref_mask_SUVR_rrrSUV_UCBJ_' subj '_Mask_VOI_SO_AtlasspaceSPM_UCBJ_' subj '_60min_90min_refVOI.nii']);
    seg_path=fullfile(FSdir,subj,'mri','aparc+aseg.nii');
    % Define thresholds, scantimes, 
    tracer='UCBJ';
%     GM_threshold       = 0.3;
%     WM_threshold       = 0.2;
    SUVR_start_time    = 60; % in min
    SUVR_end_time      = 90; % in min
    
%     % Coregister FLAIR and WMH to T1
%     [~,WMHpath]=LTNP_spm12_coregister(T1path,FLAIRpath,out_folder_FS,WMHpath);
%     
%     % Coregister WM, ref_VOI, and autoWMH to T1
%     [~,rother_image]=LTNP_spm12_coregister_reslice(T1path,T1path_cat12,out_folder_FS,{WMpath;ref_VOI_patient;autoWMHpath});
%     %rWMHpath=fullfile(out_folder_flair,['r' subj '_lesions.nii']); 
%     WMpath=rother_image{1};
%     ref_VOI_patient=rother_image{2};
%     autoWMHpath=rother_image{3};
    
    % Coregister SUV to T1
    rrrSUVnii=LTNP_spm12_coregister_reslice(T1path, rSUVnii_path, out_folder_FS); % first rigid coreg with spm

    % Add script_dir to matlab path
    addpath(rbv_script_dir)
    cd(rbv_script_dir)
    rbv_path=py.PVC_RBV_calculation_Nathalie.rbv_pvc(subj,rrrSUVnii,seg_path,out_folder_FS);
    rbv_path=char(rbv_path); 

    % Read data
%     [SUV, Vref]=LCN12_read_image(rrrSUVnii);
%     %GMimg=LCN12_read_image(GMpath,Vref);
%     WMimg=LCN12_read_image(WMpath,Vref);
%     refVOIimg=LCN12_read_image(ref_VOI_patient,Vref);
%     brain_mask=LCN12_read_image(BRAINmask_path,Vref); 
%     WMHimg=LCN12_read_image(WMHpath,Vref);
%     
%     % Threshold
%     refVOI_thresholded=refVOIimg>0.5;
%     brain_mask_thresholded=brain_mask>0.5;
%     rWMHimg_thresholded=WMHimg>0.5;
%       
%     % Use the LTNP_calculate_SUVR functions   
%     [SUVR_img, ref_mask_SO_without_manWMH, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_ucbj(SUV,refVOI_thresholded,WMimg,WM_threshold,brain_mask_thresholded,WMHimg);
%     [~, ref_mask_SO, ~,~]=LTNP_calculate_SUVR_ucbj(SUV,refVOI_thresholded,WMimg,WM_threshold,brain_mask_thresholded,'');


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



    % Save SUVR, ref mask, rWMHimg_thresholded
%     [~, SUV_name, ~]=fileparts(PETpath);
%     ref_mask_SO_path=fullfile(out_folder_suvr,['no_WMLcorrection_ref_mask_SUVR_' SUV_name '_' ref_VOI_name '_' tracer '_' subj '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min_refVOI.nii']);
%     ref_mask_SO_without_manWMH_path = fullfile(out_folder_suvr,['manualWML_masked_' 'ref_mask_SUVR_' SUV_name '_' ref_VOI_name '_' tracer '_' subj '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min_refVOI.nii']);
%     SUVR_path = fullfile(out_folder_suvr,['SUVR_' SUV_name '_WML_masked_' ref_VOI_name '_' tracer '_' subj '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min.nii']); 
%     rWMHimg_mask_path = fullfile(out_folder_FS,['rfs' subj '_lesions_binary.nii']);
%     LCN12_write_image(SUVR_img,SUVR_path,'SUVR',Vref.dt(1),Vref); 
%     LCN12_write_image(ref_mask_SO_without_manWMH,ref_mask_SO_without_manWMH_path,'SUVR_refVOI',Vref.dt(1),Vref);  
%     LCN12_write_image(ref_mask_SO,ref_mask_SO_path,'SUVR_refVOI',Vref.dt(1),Vref);  
%     LCN12_write_image(rWMHimg_thresholded,rWMHimg_mask_path,'SUVR_refVOI',Vref.dt(1),Vref);  
    
%     % Calculate WML uptake
%     VOIdetails_path='';
%     [table_autoWMH,~,~]=LTNP_VOI_stats_v8(PETpath,autoWMHpath,VOIdetails_path);
%     [table_manWMH,~,~]=LTNP_VOI_stats_v8(PETpath,rWMHimg_mask_path,VOIdetails_path);
%     [table_SO,~,~]=LTNP_VOI_stats_v8(PETpath,ref_mask_SO_path,VOIdetails_path);
%     [table_SO_without_autoWMH,~,~]=LTNP_VOI_stats_v8(PETpath,ref_mask_SO_without_autoWMH_path,VOIdetails_path);
%     [table_SO_without_manWMH,~,~]=LTNP_VOI_stats_v8(PETpath,ref_mask_SO_without_manWMH_path,VOIdetails_path);
%     
%     % Save refVOI parameters
%     if isequal(s,1)
%         table_refVOI=cell(1+length(subjects),7);
%         table_refVOI{1,1} = 'subjects'; % Create column headers for table  
%         table_refVOI{1,2} = 'mean_uptake_SO'; % Create column headers for table
%         table_refVOI{1,3} = 'nr_voxels_SO'; % Create column headers for table
%         table_refVOI{1,4} = 'mean_uptake_SO_autoWMH'; % Create column headers for table  
%         table_refVOI{1,5} = 'nr_voxels_SO_autoWMH'; % Create column headers for table
%         table_refVOI{1,6} = 'mean_uptake_SO_manWMH'; % Create column headers for table
%         table_refVOI{1,7} = 'nr_voxels_SO_manWMH'; % Create column headers for table
%         table_refVOI{1+s,1}=subj;
%         table_refVOI{1+s,2}=table_SO{2,2};
%         table_refVOI{1+s,3}=table_SO{2,7};
%         table_refVOI{1+s,4}=table_SO_without_autoWMH{2,2};
%         table_refVOI{1+s,5}=table_SO_without_autoWMH{2,7};
%         table_refVOI{1+s,6}=table_SO_without_manWMH{2,2};
%         table_refVOI{1+s,7}=table_SO_without_manWMH{2,7};
%     else
%         table_refVOI{1+s,1}=subj;
%         table_refVOI{1+s,2}=table_SO{2,2};
%         table_refVOI{1+s,3}=table_SO{2,7};
%         table_refVOI{1+s,4}=table_SO_without_autoWMH{2,2};
%         table_refVOI{1+s,5}=table_SO_without_autoWMH{2,7};
%         table_refVOI{1+s,6}=table_SO_without_manWMH{2,2};
%         table_refVOI{1+s,7}=table_SO_without_manWMH{2,7};
%     end
%     
%     % Save WMH parameters
%     if isequal(s,1)
%         table_WMH=cell(1+length(subjects),5);
%         table_WMH{1,1} = 'subjects'; % Create column headers for table  
%         table_WMH{1,2} = 'mean_uptake_autoWMH'; % Create column headers for table
%         table_WMH{1,3} = 'nr_voxels_autoWMH'; % Create column headers for table
%         table_WMH{1,4} = 'mean_uptake_manWMH'; % Create column headers for table
%         table_WMH{1,5} = 'nr_voxels_manWMH'; % Create column headers for table
%         table_WMH{1+s,1}=subj;
%         table_WMH{1+s,2}=table_autoWMH{2,2};
%         table_WMH{1+s,3}=table_autoWMH{2,7};
%         table_WMH{1+s,4}=table_manWMH{2,2};
%         table_WMH{1+s,5}=table_manWMH{2,7};
%     else
%         table_WMH{1+s,1}=subj;
%         table_WMH{1+s,2}=table_autoWMH{2,2};
%         table_WMH{1+s,3}=table_autoWMH{2,7};
%         table_WMH{1+s,4}=table_manWMH{2,2};
%         table_WMH{1+s,5}=table_manWMH{2,7};
%     end
    
end

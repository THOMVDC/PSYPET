% Define directories
tracer='UCBJ';
%WMHdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/WML/';
T1PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET2/UCBJ/';
%FLAIRdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/FLAIR/';
%PSYPETFSdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPETFASTSURFER/';
FSdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/FASTSURFER/no_biais_correction/fastsurfer_output_1_psypet1.0';
rbv_script_dir='/Volumes/LaCie/Thomas/Projects/SCRIPTS/PSYPET';
VOIdetails_path='/KUL_apps/freesurfer/FreeSurferColorLUT_VOIdetails.csv';
outdir='/Volumes/LaCie/Thomas/Projects/UCBJ_paper/UCBJ_paper_fastsurfer_pvc_brainstem';
name_overview=fullfile(outdir,[tracer '_FS_overview.xlsx']); 
name_volumes=fullfile(outdir,[tracer '_FS_volumes.xlsx']);
name_REFVOI_table=fullfile(outdir,[tracer '_FS_refVOI.xlsx']);

%subjects=dir([T1PETdir 'B0*']);
subjects={'B069','B072','B074','B076'};

for s=1:length(subjects)
    
    % Grab subject code
    %subj=subjects(s).name;
    subj=subjects{s};
    
    % Define and make output folder
    in_folder=dir(fullfile(T1PETdir,subj,'psypet_v2.0_processed_CAT12.7_*2022'));
    mkdir(fullfile(in_folder.folder,in_folder.name),'FS');
    out_folder_FS=fullfile(in_folder.folder,in_folder.name,'FS');
    ref_VOI_patient=16; %fullfile(FSdir,subj,'mri',[ref_VOI_name '.nii']); % the number corresponds to the VOI in the seg_path that will serve as ref. voi
        
    % Grab corresponding images
    rSUVnii_path=fullfile(in_folder.folder,in_folder.name,'UCBJ','REALIGNED',['rSUV_UCBJ_' subj '.nii']);
    T1path=fullfile(FSdir,subj,'mri','T1.nii');
    seg_path=fullfile(FSdir,subj,'mri','aparc+aseg.nii');
        
    % Coregister SUV to T1
    rrrSUVnii=LTNP_spm12_coregister_reslice(T1path, rSUVnii_path, out_folder_FS); % first rigid coreg with spm

    % Add script_dir to matlab path
    addpath(rbv_script_dir)
    cd(rbv_script_dir)
    rbv_path=py.PVC_RBV_calculation_Nathalie.rbv_pvc(subj,rrrSUVnii,seg_path,out_folder_FS);
    rbv_path=char(rbv_path); 

    % Read data
    [SUV, Vref]=LCN12_read_image(rbv_path);
    %GMimg=LCN12_read_image(GMpath,Vref);
    ref_mask=LCN12_read_image(seg_path,Vref);

    % Use the LTNP_calculate_SUVR functions   
    ref_value= nanmean(SUV(ref_mask==16));
    ref_volume= sum(ref_mask(:)==16);

    % Calculate SUVR
    SUVR             = SUV./ref_value;

    % Save SUVR, ref mask, rWMHimg_thresholded
    [~, SUV_name, ~]=fileparts(rbv_path);
    SUVR_path = fullfile(out_folder_FS,['SUVR_' SUV_name '_' num2str(ref_VOI_patient) '_' subj '.nii']); 
    LCN12_write_image(SUVR,SUVR_path,'SUVR',Vref.dt(1),Vref); 
    
    % Calculate stats
    [SUVR_table,~,~]=LTNP_VOI_stats_v8(SUVR_path,seg_path,VOIdetails_path);

    % Save table
    SUVR_table_path=fullfile(outdir,[tracer '_refVOI_stats.xlsx']);
    writecell(SUVR_table,SUVR_table_path,'Sheet',subj);
    
    % Add to overview
    SUVR_table(1,2)={subj};
    overview=cell(1+length(subjects),96);
    if s==1
        overview=[SUVR_table(:,1)';SUVR_table(:,2)'];
    else
        overview=[overview; SUVR_table(:,2)'];
    end
    
    % Add to volume overview 
    overview=cell(1+length(subjects),96);
    if s==1
        volumes=[SUVR_table(:,1)';SUVR_table(:,7)'];
    else
        volumes=[volumes; SUVR_table(:,7)'];
    end
    
    if s==1
        REFVOI_table=cell(1+s,3);
    end
    REFVOI_table{1+s,1}=subj;
    REFVOI_table{1+s,2}=ref_value;
    REFVOI_table{1+s,3}=ref_volume;

end

% Write overview
writecell(overview,name_overview)
writecell(volumes,name_volumes)
writecell(REFVOI_table,name_REFVOI_table)



% better to take brain.nii as T1 than nu.nii for bringing the SO to subject
% <-> but quid for corgestritaiton pet/t1 : is nu.nii better or brain.nii ?
% Doesnot look like so; but T1.nii did well. T1 is also the last being
% processed before aparc+asegDKT, so we'll take T1.nii
% aparc.DKTatlas+aseg.deep.withCC.mgz

%ect='post';

% Define input and output dir
%PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/pet_output/';
%PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/psypet_fs_output';
%PETdir='/Volumes/LaCie/Margot/aline/UCBJ_psypet1';
%PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/MK62_4s';
PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/UCBJ_4s';
%PETdir=['/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/UCBJ_4s_ECT/' ect];
%PETdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_cat12_output_SOnoWML_as_rr';
%T1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/fastsurfer_output';
%T1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/T1/DICOM';
T1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/input';
%T1dir=['/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4_ECT/fs/imaging/'  ect];
%T1dir_atlas='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/cat12_output/';
T1dir_atlas='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/psycat';
%T1dir='/Volumes/LaCie/Margot/aline/T1_fastsurfer/fastsurfer_output';
%deffolder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/cat12_output/mri';
%outfolder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/pet_output';
%outfolder=['/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4_ECT/fs/imaging/'  ect];
outfolder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_cat12_output_SO_mci';
%outfolder_xls=['/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4_ECT/fs/stats/'  ect];
outfolder_xls='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_cat12_stats_SO_mci';
refVOI='/Volumes/LaCie/Thomas/Projects/SCRIPTS/PSYPET/templates/TemplateSO_78subj_SPM_Juli2020.nii';
%atlas='neuromorphometrics';
%refVOI=[13,14,37,38,39,40,41,42];

% Loop through subjects B066,B070, B073, B077
subjects=dir(fullfile(PETdir, 'MCI*'));
%subjects={'B066','B070','B073','B077'};

%% Launcher for CAT12
for s=1:length(subjects)
    
    % Grab subject code
    subj=subjects(s).name;
    %subj=subjects{s};
    
    % Grab input folder
    %T1_subj=dir(fullfile(T1dir,subj,'psypet_v2.0_processed_CAT12.7_*2022','ANAT'));         
    %PET=fullfile(PETdir,subj,['rSUV_' subj '.nii']);
    PET=fullfile(PETdir,subj);
    %T1_mgz=fullfile(T1dir,subj,'mri','T1.mgz'); % t1 fastsurfer;
    %T1=fullfile(T1_subj(1).folder,T1_subj(1).name,['accT1_' subj '.nii']);
    %T1=fullfile(T1dir,['accT1_' subj '.nii']);
    T1=fullfile(T1dir,subj);
    %atlas_mgz=fullfile(T1dir_atlas,subj,'mri','aparc.DKTatlas+aseg.deep.withCC.mgz'); % aparc aseg fastsurfer
    %atlas=fullfile(T1dir_atlas,subj,['rbv_segm_neuromorphometrics_accT1_' subj '.nii']);
    %atlas=fullfile(T1dir_atlas,subj,['rbv_segm_skull_neuromorphometrics_accT1_' subj '.nii']);
    %atlas='neuromorphometrics';
    outfolder_subj=fullfile(outfolder,subj);
    cd(outfolder);
    mkdir(subj);
    
    % Convert mgz to nii
    %[T1]=LTNP_mgz2nii(T1_mgz,outfolder_subj);
    %[atlas]=LTNP_mgz2nii(atlas_mgz,outfolder_subj);
    
    % Calculate voxelsize
    %voxelsize=LTNP_get_voxelsize(T1);
    %voxelsize=1;

    % Make rr ready
    %[~,invdef]=LTNP_calc_deformation_field(T1,outfolder_subj);
    %invdef=fullfile(deffolder,['iy_accT1_' subj '.nii']);
    %[refVOI]=LTNP_spm12_warp_ROI(rr,invdef,outfolder_subj,voxelsize);

    % Run psypet4
    [SUVR_path, SUVR_table_path, SUV_rr_table_path]=psypet(subj, T1, PET, refVOI, atlas, outfolder_subj);
    
end


% All excels in one for cat12
VOIdetails='/Users/mlaroy0/spm12/toolbox/neuromorphometrics_144.csv';
%VOIdetails='/Users/mlaroy0/spm12/toolbox/neuromorphometrics_145.csv';
VOIdet=readcell(VOIdetails);
SUVR_table_all=fullfile(outfolder_xls,'SUVR_noPVC_stats.xlsx');
nr_vois=length(VOIdet);
for s=1:length(subjects)
    subj=subjects(s).name;
    %SUVR_table_path=fullfile(outfolder,subj,'SUVR_SUV_PET_PVC_RBV_65mm_in_seg_wTemplateSO_78subj_SPM_Juli2020.xlsx');
    SUVR_table_path=fullfile(PETdir,subj,['SUVR_rrrSUV_' subj '_refVOI_without_WML.xlsx']);
    SUVR_table=readcell(SUVR_table_path);
    SUVR_table_update=cell(nr_vois+1,size(SUVR_table,2));
    SUVR_table(1,1)={subj};
    SUVR_table_update(1,:)=SUVR_table(1,:);
    for v = 1:nr_vois
        ID=num2str(VOIdet{v,1}); % grabs the ROIids
        SUVR_table_update(v+1,1)=VOIdet(v,3);
        if ismember(ID,SUVR_table(:,1)) % only consider values from the atlas image
            [~,index]=ismember(ID,SUVR_table(1:end,1));
            SUVR_table_update(v+1,2:end)=SUVR_table(index,2:end);
        end
    end
    writecell(SUVR_table_update,SUVR_table_all,'sheet',subj)
end

%% launcher for freesurfer of fastsurfer
for s=1:length(subjects)
    
    % Grab subject code
    subj=subjects(s).name;
    %subj=subjects{s};
    
    % Grab input folder
    %T1_subj=dir(fullfile(T1dir,subj,'psypet_v2.0_processed_CAT12.7_*2022','ANAT'));         
    %PET=fullfile(PETdir,subj,['rSUV_' subj '.nii']);
    PET=fullfile(PETdir,subj);
    T1=fullfile(T1dir,subj,'T1.nii'); % t1 fastsurfer;
    %T1=fullfile(T1_subj(1).folder,T1_subj(1).name,['accT1_' subj '.nii']);
    %T1=fullfile(T1dir,['accT1_' subj '.nii']);
    %T1=fullfile(T1dir,subj);
    atlas=fullfile(T1dir,subj,'aparc+aseg.nii'); % aparc aseg fastsurfer
    %atlas=fullfile(T1dir_atlas,subj,['rbv_segm_neuromorphometrics_accT1_' subj '.nii']);
    %atlas=fullfile(T1dir_atlas,subj,['rbv_segm_skull_neuromorphometrics_accT1_' subj '.nii']);
    %atlas='neuromorphometrics';
    outfolder_subj=fullfile(outfolder,subj);
    %cd(outfolder);
    %mkdir(subj);
    
    % Convert mgz to nii
    %[T1]=LTNP_mgz2nii(T1_mgz,outfolder_subj);
    %[atlas]=LTNP_mgz2nii(atlas_mgz,outfolder_subj);
    
    % Calculate voxelsize
    voxelsize=LTNP_get_voxelsize(T1);
    %voxelsize=1;

    % Make rr ready
    [~,invdef]=LTNP_calc_deformation_field(T1,outfolder_subj);
    %invdef=fullfile(deffolder,['iy_accT1_' subj '.nii']);
    [refVOI_warped]=LTNP_spm12_warp_ROI(refVOI,invdef,outfolder_subj,voxelsize);

    % Run psypet4
    [SUVR_path, SUVR_table_path, SUV_rr_table_path]=psypet(subj, T1, PET, refVOI_warped, atlas, outfolder_subj);
    
end

% All excels in one for fs
%VOIdetails='/KUL_apps/freesurfer/FastSurferColorLUT_VOIdetails_100.csv'; %
%fastsurfer
VOIdetails='/KUL_apps/freesurfer/FreeSurferColorLUT.txt'; % freesurfer
VOIdet=readcell(VOIdetails);
SUVR_table_all=fullfile(outfolder_xls,['SUVR_PVC_stats_' ect '.xlsx']);
nr_vois=length(VOIdet);
for s=1:length(subjects)
    subj=subjects(s).name;
    SUVR_table_path=fullfile(outfolder,subj,'SUVR_SUV_PET_PVC_RBV_65mm_in_seg_wTemplateSO_78subj_SPM_Juli2020.xlsx');
    SUVR_table=readcell(SUVR_table_path);
    SUVR_table_update=cell(nr_vois+1,size(SUVR_table,2));
    SUVR_table(1,1)={subj};
    SUVR_table_update(1,:)=SUVR_table(1,:);
    for v = 1:nr_vois
        ID=num2str(VOIdet{v,1}); % grabs the ROIids
        SUVR_table_update(v+1,1)=VOIdet(v,2); % grab the name corresponding to the ROIid
        if ismember(ID,SUVR_table(:,1)) % only consider values from the atlas image
            [~,index]=ismember(ID,SUVR_table(1:end,1));
            SUVR_table_update(v+1,2:end)=SUVR_table(index,2:end);
        end
    end
    row_is_empty = any(cellfun(@isempty, SUVR_table_update), 2); % find rows with empty values
    SUVR_table_update(row_is_empty,:) = []; % delete these rows
    writecell(SUVR_table_update,SUVR_table_all,'sheet',subj)
end



    
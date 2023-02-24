%% CAT12 SUVR stats
atlas_name='neuromorphometrics';
VOIdetails_path=fullfile('/Users/mlaroy0/spm12/toolbox/cat12/templates_volumes', [atlas_name '.csv']);
tracer='UCBJ';
%petdir=fullfile('/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/', [tracer '_4s']); % directory with processed subject folders
petdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/UCBJ_4s';
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/T1';
outdir='/Volumes/LaCie/Thomas/Projects/UCBJ_paper/UCBJ_paper_20lld_2_pvc';
name_excel=fullfile(outdir,[tracer '_PVC_stats.xlsx']); % give the desired name and path to the excel output
cd(petdir)
l=dir('*0*');

% Choose in which matter you look for
GM=1;
WM=0;
CSF=0;
WMHC=0;

for i=1:length(l)
    
    % Grab subject
    subj=l(i).name;
    disp(subj);
    
    % Grab necessary subject files
    Tname=['accT1_' subj '.nii'];

    GMmask_path=fullfile(t1dir,subj,'MASKS',['mask_p1' Tname]); % change path construction to grab the desired GM mask
    WMmask_path=fullfile(t1dir,subj,'MASKS',['mask_p2' Tname]);
    CSFmask_path=fullfile(t1dir,subj,'MASKS',['mask_p3' Tname]);
    WMHmask_path=fullfile(t1dir,subj,'MASKS',['mask_p7' Tname]);
%   GMmask_path=fullfile(t1dir,subj,'MASKS',['mask_mwp1' Tname]); % change path construction to grab the desired GM mask
%   WMmask_path=fullfile(t1dir,subj,'MASKS',['mask_mwp2' Tname]);
%   CSFmask_path=fullfile(t1dir,subj,'MASKS',['mask_mwp3' Tname]);
%   WMHmask_path=fullfile(t1dir,subj,'MASKS',['mask_mwp7' Tname]);
    
    % If you don't want PVC
    %SUVR_dir=dir(fullfile(petdir,subj,'psypet*',tracer,'SUVR',['SUVR_rrrSUV_*' subj '*.nii'])); % change path construction to grab the desired SUVR
    %atlas_path=fullfile(t1dir,subj,'CAT12','mri_atlas',[atlas_name '_' Tname]);

    % If you want RBV
    atlas_path=fullfile(t1dir,subj,'RBV_ATLASSES',['rbv_segm_' atlas_name '_' Tname]); 
    SUVR_dir=dir(fullfile(petdir,subj,'psypet*',tracer,'SUVR',['SUVR_PVC_RBV2_' subj '*.nii']));
    
    % Check if there is only one PET file
    if length(SUVR_dir)==1
        SUVR_path=fullfile(SUVR_dir(1).folder,SUVR_dir(1).name);
    else
        fprintf(subj,'0 or more than one SUVR image available \n');
    end
    
    % Calculate stats
    if GM==1
        [SUVR_GM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,GMmask_path);
    else
        SUVR_GM_table=[];
    end
    if WM==1
        [SUVR_WM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,WMmask_path);
    else
        SUVR_WM_table=[];
    end
    if CSF==1
        [SUVR_CSF_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,CSFmask_path);
    else
        SUVR_CSF_table=[];
    end
    if WMHC==1
        [SUVR_WMH_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,WMHmask_path);
    else 
        SUVR_WMH_table=[];
    end
    % Concatenate stats
    SUVR_table=horzcat(SUVR_GM_table(:,1:end),SUVR_WM_table(:,2:end),SUVR_CSF_table(:,2:end),SUVR_WMH_table(:,2:end));
    
    % Write table
    SUVR_table{1,1}=subj;
    writecell(SUVR_table,name_excel,'sheet',subj)
    
end

%% AMY positivity stats%;
tracer='UCBJ';
%petdir=['/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/' tracer '_4s'];
petdir='/Volumes/LaCie/Thomas/Projects/RETINAL_IMAGING/DATA/CENTILOID/';
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/T1';
%name_excel='/Volumes/LaCie/Thomas/Projects/UCBJ_paper/T1_compVOI.xlsx';
name_excel='/Volumes/LaCie/Thomas/Projects/RETINAL_IMAGING/DATA/CENTILOID/centiloid_compVOI_values.xlsx';
'B046/processed/clSUVR_cerebellum_90min_120min.nii
%name_excel='/Users/tvdcas1/Desktop/rbv_amy_pos.xls';
atlas_path='/Volumes/LaCie/Thomas/Projects/SCRIPTS/AMYWM/LCN_Dupont/resl_to_wmc1_SPM12_aal_composite_cortical.img'; % mask or atlas
VOIdetails_path='/Volumes/LaCie/Thomas/Projects/SCRIPTS/AMYWM/LCN_Dupont/resl_to_wmc1_SPM12_aal_composite_cortical.csv'; % csv with voi details of the given atlas/mask
subjects=dir(fullfile(petdir, '*'));
for i=1:length(subjects)
    subj=subjects(i).name;
    disp(subj)
    SUVR_path=fullfile(petdir,subj,'processed','clSUVR_cerebellum_90min_120min.nii');
    %SUVR_dir=dir(fullfile(petdir,subj,'psypet*',tracer,'SUVR',['SUVR_wrrrSUV_*' subj '*.nii'])); % change path construction to grab the desired SUVR
    %SUVR_dir=dir(fullfile(petdir,subj,'psypet*',tracer,'SUVR',['SUVR_wPVC_RBV_' subj '*.nii'])); % rbv one
    %SUVR_path=fullfile(SUVR_dir(1).folder,SUVR_dir(1).name);
    GMmask_path=fullfile(t1dir,subj,'CAT12','mri',['mwp1accT1_' subj '.nii']); % change path construction to grab the desired GM mask
    %SUVR_path=GMmask_path;
    [SUVR_GM_table,~,~]=LTNP_VOI_stats(SUVR_path,atlas_path,VOIdetails_path,GMmask_path);
    if i==1
        SUVR_GM_table{1,1}=SUVR_GM_table{2,1};
        SUVR_GM_table{2,1}=subj;
        writecell(SUVR_GM_table,name_excel)
    else
        SUVR_GM_table{2,1}=subj;
        writecell(SUVR_GM_table(2,:),name_excel,'WriteMode','append');
    end
end

%% ANATOMICAL CAT12 stats
atlas_name='neuromorphometrics';
VOIdetails_path=fullfile('/Users/mlaroy0/spm12/toolbox/cat12/templates_volumes', [atlas_name '.csv']);
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/T1';
name_excel='/Volumes/LaCie/Thomas/Projects/UCBJ_paper/UCBJ_paper_20lld/T1_stats.xlsx';
cd(t1dir)
l=dir('*0*'); 
for i=1:length(l)
    
    % grab subject
    subj=l(i).name;
    Tname=['accT1_' subj '.nii'];
    disp(subj);
    
    % grab segmentations
    GM_path=fullfile(t1dir,subj,'CAT12','mri',['p1' Tname]); % change path construction to grab the desired GM mask
    WM_path=fullfile(t1dir,subj,'CAT12','mri',['p2' Tname]);
    CSF_path=fullfile(t1dir,subj,'CAT12','mri',['p3' Tname]);
    WMH_path=fullfile(t1dir,subj,'CAT12','mri',['p7' Tname]);
    
    % Grab thickness map
    thickness_path=fullfile(t1dir,subj,'THICKNESS',['thickness_p1' Tname]); 
    
    % grab masks
    GMmask_path=fullfile(t1dir,subj,'MASKS',['mask_p1' Tname]); % change path construction to grab the desired GM mask
    WMmask_path=fullfile(t1dir,subj,'MASKS',['mask_p2' Tname]);
    CSFmask_path=fullfile(t1dir,subj,'MASKS',['mask_p3' Tname]);
    WMHmask_path=fullfile(t1dir,subj,'MASKS',['mask_p7' Tname]);
    
    % grab atlas
    atlas_path=fullfile(t1dir,subj,'CAT12','mri_atlas',[atlas_name '_' Tname]);
    
    % calculate stats
    [GM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,GM_path,GMmask_path);
    [WM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,WM_path,WMmask_path);
    [CSF_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,CSF_path,CSFmask_path);
    [WMH_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,WMH_path,WMHmask_path);
    [thickness_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,thickness_path);
    
    % Concatenate tables
    stats_p1p2p3p7_thick=horzcat(GM_table(:,1:end),WM_table(:,2:end),CSF_table(:,2:end),WMH_table(:,2:end),thickness_table(:,2:end));
    
    % Write to excel
    stats_p1p2p3p7_thick{1,1}=subj;
    writecell(stats_p1p2p3p7_thick,name_excel,'sheet',subj)

end

%% Meningeal stats
for i=1:length(l)
    
    % Grab subject
    subj=l(i).name;
    disp(subj)
    
    % Grab necessary subject files
    Tname=['accT1_' subj '.nii'];
    rbv_atlas_path=fullfile(t1dir,subj,'RBV_ATLASSES',['rbv_segm_' atlas_name '_' Tname]); 
    VOIdetails_rbv=fullfile(t1dir,subj,'RBV_ATLASSES','VOIdetails_rbv.csv');
    SUVR_dir=dir(fullfile(petdir,subj,'psypet*','UCBJ','SUVR',['SUVR_rrrSUV_*' subj '*.nii'])); % change path construction to grab the desired SUVR
    
    % Check if there is only one PET file
    if length(SUVR_dir)==1
        SUVR_path=fullfile(SUVR_dir(1).folder,SUVR_dir(1).name);
    else
        fprintf(subj,'\n');
    end
    
    % Calculate stats
    [SUVR_MENG_table,~,~]=LTNP_VOI_stats_v7(rbv_atlas_path,VOIdetails_rbv,SUVR_path);
        
    % Write table
    if i==1
        SUVR_table{1,1}=SUVR_table{2,1};
        SUVR_table{2,1}=subj;
        writecell(SUVR_table,name_excel)
    else
        SUVR_table{2,1}=subj;
        writecell(SUVR_table(2,:),name_excel,'sheet',subj);
    end
end

%% CAT12 tissue stats
tracer='UCBJ';
petdir=fullfile('/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/', [tracer '_4s']); % directory with processed subject folders
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/T1';
outdir='/Volumes/LaCie/Thomas/Projects/UCBJ_paper/UCBJ';
name_excel=fullfile(outdir,[tracer '_4s_wholebrain.xlsx']); % give the desired name and path to the excel output
cd(petdir)
l=dir('*0*');

% Choose in which matter you look for
tissue='GM'; %WM or CSF or WMHC 


for i=1:length(l)
    
    % Grab subject
    subj=l(i).name;
    disp(subj);
    
    % Grab necessary subject files
    Tname=['accT1_' subj '.nii'];

    SUVR_dir=dir(fullfile(petdir,subj,'psypet*',tracer,'SUVR',['SUVR_rrrSUV_*' subj '*.nii'])); % change path construction to grab the desired SUVR
    %SUVR_dir=dir(fullfile(petdir,subj,'psypet*',tracer,'SUVR',['SUVR_PVC_RBV_' subj '*.nii']));
    
    % Check if there is only one PET file
    if length(SUVR_dir)==1
        SUVR_path=fullfile(SUVR_dir(1).folder,SUVR_dir(1).name);
    else
        fprintf(subj,'0 or more than one SUVR image available \n');
    end
    
    % Calculate stats
    if strcmp(tissue,'GM')
        atlas_path=fullfile(t1dir,subj,'MASKS',['mask_p1accT1_' subj '.nii']);
        ROIid = 1;
        ROIabbr = {[tissue 'tot']};
        ROIname = {[tissue 'total']};
        VOIdetails_table=table(ROIid,ROIabbr,ROIname);
        VOIdetails_path=fullfile(outdir,[tissue 'tot.csv']);
        writetable(VOIdetails_table,VOIdetails_path);
        [SUVR_GM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
        SUVR_GM_table{2,1}=subj;
        if i==1
            SUVR_table=SUVR_GM_table;
        else
            SUVR_table=vertcat(SUVR_table, SUVR_GM_table(2,1:end));
        end
    end
    if strcmp(tissue,'WM')
        atlas_path=fullfile(t1dir,subj,'MASKS',['mask_p2accT1_' subj '.nii']);
        ROIid = 1;
        ROIabbr = {[tissue 'tot']};
        ROIname = {[tissue 'total']};
        VOIdetails_table=table(ROIid,ROIabbr,ROIname);
        VOIdetails_path=fullfile(outdir,[tissue 'tot.csv']);
        writetable(VOIdetails_table,VOIdetails_path);
        [SUVR_WM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
        SUVR_WM_table{2,1}=subj;
        if i==1
            SUVR_table=SUVR_WM_table;
        else
            SUVR_table=vertcat(SUVR_table, SUVR_WM_table(2,1:end));
        end
    end
    if strcmp(tissue,'CSF')
        atlas_path=fullfile(t1dir,subj,'MASKS',['mask_p3accT1_' subj '.nii']);
        ROIid = 1;
        ROIabbr = {[tissue 'tot']};
        ROIname = {[tissue 'total']};
        VOIdetails_table=table(ROIid,ROIabbr,ROIname);
        VOIdetails_path=fullfile(outdir,[tissue 'tot.csv']);
        writetable(VOIdetails_table,VOIdetails_path);
        [SUVR_CSF_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
        SUVR_CSF_table{2,1}=subj;
        if i==1
            SUVR_table=SUVR_CSF_table;
        else
            SUVR_table=vertcat(SUVR_table, SUVR_CSF_table(2,1:end));
        end
    end
    if strcmp(tissue,'WMHC')
        atlas_path=fullfile(t1dir,subj,'MASKS',['mask_p7accT1_' subj '.nii']);
        ROIid = 1;
        ROIabbr = {[tissue 'tot']};
        ROIname = {[tissue 'total']};
        VOIdetails_table=table(ROIid,ROIabbr,ROIname);
        VOIdetails_path=fullfile(outdir,[tissue 'tot.csv']);
        writetable(VOIdetails_table,VOIdetails_path);
        [SUVR_WMH_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
        SUVR_WMH_table{2,1}=subj;
        if i==1
            SUVR_table=SUVR_WMH_table;
        else
            SUVR_table=vertcat(SUVR_table, SUVR_WMH_table(2,1:end));
        end
    end
end

% Write table
writetable(SUVR_table,name_excel)

%% Fastsurfer SUVR stats
atlas_name='fastsurfer';
VOIdetails_path='/KUL_apps/freesurfer/FreeSurferColorLUT_VOIdetails_100.csv'; %fullfile('/Users/mlaroy0/spm12/toolbox/cat12/templates_volumes', [atlas_name '.csv']);
%VOIdetails_path='';
tracer='UCBJ';
petdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_fs_output/';
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_fs_output/';
outdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_fs_output/';
name_excel=fullfile(outdir,[tracer '_FS_stats.xlsx']); % give the desired name and path to the excel output
name_overview=fullfile(outdir,[tracer '_FS_overview.xlsx']); 
name_volumes=fullfile(outdir,[tracer '_FS_volumes.xlsx']);
cd(petdir)
l=dir(fullfile(petdir,'B*'));

% Choose in which matter you look for
for i=1:length(l)
    
    
    % Grab subject
    subj=l(i).name;
    disp(subj);
    
    % Name excel
    name_excel=fullfile(petdir,subj,['SUVR_rrrSUV_' subj '.xlsx']);


    % Grab files
    atlas_path=fullfile(t1dir,subj,'aparc.DKTatlas+aseg.deep.withCC.nii'); 
    %SUVR_path=fullfile(petdir,['pvc_' subj '_PET_PVC_RBV_65mm_in_seg.nii']);
    SUVR_path=fullfile(petdir,subj,['SUVR_rrrSUV_' subj '.nii']);
    
    % Calculate stats
    %[SUVR_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
    [SUVR_table,~,~]=LTNP_VOI_stats(SUVR_path,atlas_path,VOIdetails_path);
    
    % Write table
    writecell(SUVR_table,name_excel) %,'sheet',subj)
    
end

    
    % Add to overview
    SUVR_table(1,2)={subj};
    if i==1
        overview=[SUVR_table(:,1)';SUVR_table(:,2)'];
    else
        overview=[overview; SUVR_table(:,2)'];
    end
    
    % Add to volume overview 
    SUVR_table(1,2)={subj};
    if i==1
        volumes=[SUVR_table(:,1)';SUVR_table(:,7)'];
    else
        volumes=[volumes; SUVR_table(:,7)'];
    end
end

% Write overview
writecell(overview,name_overview)
writecell(volumes,name_volumes)

% for i=1:length(l)
%     subj=l(i).name(20:23);
%     T=readtable(name_excel,'sheet',subj);
%     if i==1
%         volumes=[rows2vars(T(:,1));rows2vars(T(:,7))];
%     else
%         volumes=[volumes; rows2vars(T(:,7))];
%     end
%     volumes(1+i,1)={subj};
% end

%% ANATOMICAL Fastsurfer stats -> faster with fs_mr_stats2xls.py and fs_mr_xls2plot.py
VOIdetails_path='/KUL_apps/freesurfer/FreeSurferColorLUT_VOIdetails.csv'; %fullfile('/Users/mlaroy0/spm12/toolbox/cat12/templates_volumes', [atlas_name '.csv']);
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/fastsurfer_output';
outdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/fastsurfer_stats';
name_excel=fullfile(outdir,'FS_stats.xlsx'); % give the desired name and path to the excel output
l=dir(fullfile(t1dir,'*0*'));

% Choose in which matter you look for
volumes=cell(length(l)+1,96);
for i=1:length(l)
    
    % Grab subject
    subj=l(i).name;
    disp(subj);

    % Grab files
    atlas_path=fullfile(t1dir,subj,'mri','aparc.DKTatlas+aseg.deep.withCC.nii'); 
    %SUVR_path=fullfile(petdir,['pvc_' subj '_PET_PVC_RBV_65mm_in_seg.nii']);
    %SUVR_path=fullfile(l(i).folder,l(i).name);
    
    % Calculate stats
    %[SUVR_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
    [T1_table,~,~]=LTNP_VOI_stats(atlas_path,atlas_path,VOIdetails_path);
    
    % Write table
    %writecell(T1_table,name_excel,'sheet',subj)
       
    % Add to volume overview 
    T1_table(1,7)={subj};
    if i==1
        volumes=[T1_table(:,1)';T1_table(:,7)'];
    else
        volumes(1+i,:)=T1_table(:,7)';
    end
end

% Write overview
writecell(volumes,name_excel)

%% ANATOMICAL CAT12 stats psypet 4
VOIdetails_path='/Users/mlaroy0/spm12/toolbox/neuromorphometrics_144.csv'; %fullfile('/Users/mlaroy0/spm12/toolbox/cat12/templates_volumes', [atlas_name '.csv']);
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/psycat/';
outdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/cat12_stats';
name_excel=fullfile(outdir,'volumetry_stats3.xlsx'); % give the desired name and path to the excel output
l=dir(fullfile(t1dir,'*0*'));

volumes=cell(length(l)+1,145);
for i=1:length(l)
    
    % Grab subject
    subj=l(i).name;
    disp(subj);

    % Grab files
    atlas_path=fullfile(t1dir,subj,['rbv_segm_neuromorphometrics_accT1_' subj '.nii']); 
    %SUVR_path=fullfile(petdir,['pvc_' subj '_PET_PVC_RBV_65mm_in_seg.nii']);
    %SUVR_path=fullfile(l(i).folder,l(i).name);
    
    % Calculate stats
    %[SUVR_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path);
    [T1_table,~,~]=LTNP_VOI_stats(atlas_path,atlas_path,VOIdetails_path);
    
    % Write table
    %writecell(T1_table,name_excel,'sheet',subj)
       
    % Add to volume overview 
    T1_table(1,7)={subj};
    if i==1
        volumes=[T1_table(:,1)';T1_table(:,7)'];
    else
        volumes(1+i,:)=T1_table(:,7)';
    end
end

% Write overview
writecell(volumes,name_excel)



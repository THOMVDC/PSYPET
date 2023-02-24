
%% Settings

% Define input directories
%fdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/cat12_output/mri/';
%fdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/MK62/output/psypet_cat12_output/';
fdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_cat12_output/';
t1dir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/cat12_output/mri/';

EXCEL_DIR='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/C2_tau_flut_ucbj.xlsx';
sheet='FDP';

% Define output directories
OUTPUT_DIR='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/vba';
%OUTPUT_DIR='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/vbm';

% Load spread sheet
con=readtable(EXCEL_DIR,'Sheet','CON');
pat=readtable(EXCEL_DIR,'Sheet','LLD');

% Get subjectcodes
SUBJECTCODES_1=table2cell(con(:,1));
SUBJECTCODES_2=table2cell(pat(:,1));

% Remove empty cells from array
SUBJECTCODES_1=SUBJECTCODES_1(~cellfun('isempty',SUBJECTCODES_1));
SUBJECTCODES_2=SUBJECTCODES_2(~cellfun('isempty',SUBJECTCODES_2));
SUBJECTCODES=[SUBJECTCODES_1;SUBJECTCODES_2];

% Give a name to the mask that will be created
MASK_NAME=fullfile(OUTPUT_DIR,'GM_mask03.nii');

% Define settings
clin_data_name_1='age';
%clin_data_name_1=[]; % if you don't want to add clin_data
%clin_data_name_2='TIV';
clin_data_name_2=[]; % if you don't want to add clin_data
contrast_name='Con > LLD';
contrast_name_inv='LLD > Con';
convec=[1 -1];
convec_inv=[-1 1]; 
sm=8; % smooth kernel
sprefix='s'; % smooth prefix
voxelsize=1.5;

% Find number of subjects
nr_subj=length(SUBJECTCODES);
nr_subj_1=length(SUBJECTCODES_1);
nr_subj_2=length(SUBJECTCODES_2);

%% Proces images

% Grab files
files=cell(nr_subj,1);
for i=1:nr_subj
    %files{i} = fullfile(fdir,['mwp1accT1_' SUBJECTCODES{i} '.nii']);
    fname='SUVR_SUV_PET_PVC_RBV_65mm_in_seg.nii';
    files{i}=fullfile(OUTPUT_DIR,[SUBJECTCODES{i} '_' fname]);
    copyfile(fullfile(fdir,SUBJECTCODES{i},fname),files{i})
end

% warp
for i=1:nr_subj
    deformation_field=fullfile(t1dir,['y_accT1_' SUBJECTCODES{i} '.nii']);
    files{i}=LTNP_spm12_warp(files{i},deformation_field,OUTPUT_DIR,voxelsize);
end

% Smooth
matlabbatch={};
matlabbatch{1}.spm.spatial.smooth.data = files;
matlabbatch{1}.spm.spatial.smooth.fwhm = [sm sm sm]; 
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = [sprefix num2str(sm)];

spm_jobman('run',matlabbatch);

%% Create 2 lists of s8 files for VBM
files_1=cell(nr_subj_1,1);
for i = 1:nr_subj_1
    %files_1{i,1} = [fullfile(fdir,[sprefix num2str(sm) 'mwp1accT1_' SUBJECTCODES_1{i} '.nii']) ',1']; 
    files_1{i,1} = [fullfile(OUTPUT_DIR,[sprefix num2str(sm) 'w' SUBJECTCODES_1{i} '_SUVR_SUV_PET_PVC_RBV_65mm_in_seg.nii']) ',1'];
end
files_2=cell(nr_subj_2,1);
for i = 1:nr_subj_2
    %files_2{i,1} = [fullfile(fdir,[sprefix num2str(sm) 'mwp1accT1_' SUBJECTCODES_2{i} '.nii']) ',1']; 
    files_2{i,1} = [fullfile(OUTPUT_DIR,[sprefix num2str(sm) 'w' SUBJECTCODES_1{i} '_SUVR_SUV_PET_PVC_RBV_65mm_in_seg.nii']) ',1'];
end

%% Create clin_data_1
if isempty(clin_data_name_1)
    clin_data_1=[];
else
    T=readtable(EXCEL_DIR,'Sheet',sheet);
    clin_data_1=zeros((nr_subj),1);
    for j = 1:height(T)
        for i = 1:nr_subj
            if isequal(T{j,1},SUBJECTCODES(i))
                clin_data_1(i,1)=T{j,clin_data_name_1};
                disp([SUBJECTCODES(i) ' : ' T{j,clin_data_name_1}])
            end
        end
    end
end

%% Create clin_data_2
if isempty(clin_data_name_2)
    clin_data_2=[];
else
%     clin_data_2=zeros((nr_subj),1);
%     for i = 1:nr_subj_1
%         cat12_dir=dir([fdir '/' SUBJECTCODES_1{i} '/CAT12/report/cat_*.mat']);
%         cat12_tiv=load(fullfile(cat12_dir(1).folder,cat12_dir(1).name));
%         clin_data_2(i,1)=cat12_tiv.S.subjectmeasures.vol_TIV;
%     end
%     for i = 1:nr_subj_2
%         cat12_dir=dir([fdir '/' SUBJECTCODES_1{i} '/CAT12/report/cat_*.mat']);
%         cat12_tiv=load(fullfile(cat12_dir(1).folder,cat12_dir(1).name));
%         clin_data_2(nr_subj_1+i,1)=cat12_tiv.S.subjectmeasures.vol_TIV;
%     end
    T=readtable(EXCEL_DIR,'Sheet',sheet);
    clin_data_2=zeros((nr_subj),1);

    for j = 1:height(T)
        for i = 1:nr_subj
            if isequal(T{j,1},SUBJECTCODES(i))
                clin_data_2(i,1)=T{j,clin_data_name_2};
                disp([SUBJECTCODES(i) ' : ' T{j,clin_data_name_2}])
            end
        end
    end
end

%% Launch VBM
LTNP_vbm_ttest(files_1, files_2, MASK_NAME, clin_data_1, clin_data_name_1, clin_data_2, clin_data_name_2, contrast_name,contrast_name_inv, convec, convec_inv, OUTPUT_DIR)
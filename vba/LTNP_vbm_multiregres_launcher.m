%% Settings

% Define directories and subjectcodes
mdir='';
fdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/MK62/output/vba';
EXCEL_DIR='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/C2_tau_flut_ucbj.xlsx';
OUTPUT_DIR='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/MK62/output/vba/regress';
MASK_NAME='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/MK62/output/vba/ttest/GM_mask03.nii';

% Define settings
clin_data_1_name='gds_score';
contrast_name='Positive correlations';
contrast_inv_name='Negative correlations';
contrast=[0 1];
contrast_inv=[0 -1];
sm=8; % smooth kernel
sprefix='s'; % smooth prefix

% Load spread sheet
clintable=readtable(EXCEL_DIR,'Sheet','LLD');
%SUBJECTCODES=table2cell(clintable(:,1));
SUBJECTCODES={'B011';'B019';'B021';'B027';'B032';'B034';'B042';'B043';'B045';'B048';'B052';'B059';'B066';'B069';'B070';'B073';'B074';'B075';'B077';'B078';};

% Find number of subjects
nr_subj=length(SUBJECTCODES);

%% Smooth images

% Grab files
files=cell(nr_subj,1);
for i=1:nr_subj
    %files{i} = fullfile(fdir,SUBJECTCODES{i},'CAT12','mri',['mwp1accT1_' SUBJECTCODES{i} '.nii']); 
    files{i} = fullfile(fdir,'s8wSUVR',['SUVR_wrrrSUV_UCBJ_' SUBJECTCODES{i} '_Mask_VOI_SO_AtlasspaceSPM_UCBJ_' SUBJECTCODES{i} '_60min_90min.nii']);
end

% Smooth
matlabbatch={};
matlabbatch{1}.spm.spatial.smooth.data = files;
matlabbatch{1}.spm.spatial.smooth.fwhm = [sm sm sm]; 
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = [sprefix num2str(sm)];

spm_jobman('run',matlabbatch);

%% Create list of wc1T1 as input for imcalc
if isfile(MASK_NAME) 
    maskfiles=[];
else
    maskfiles=cell(nr_subj,1);
    for i = 1:nr_subj
        file=dir([mdir '/wT1/wc1T1ac_' char(SUBJECTCODES(i)) '.nii']);
        maskfiles{i,1} = [fullfile(file.folder,file.name) ',1']; 
    end
end

%% Create scan list 1
scanfiles=cell(nr_subj,1);
for i = 1:nr_subj
    %scanfiles{i,1} = [fullfile(fdir,[sprefix num2str(sm) 'SUVR_wrrrSUV_UCBJ_' SUBJECTCODES{i} '_Mask_VOI_SO_AtlasspaceSPM_UCBJ_' SUBJECTCODES{i} '_60min_90min.nii']) ',1'];
    scanfiles{i,1} = [fullfile(fdir,'s8wSUVR',[sprefix num2str(sm) 'w' SUBJECTCODES{i} '_SUVR_SUV_PET_PVC_RBV_65mm_in_seg.nii']) ',1'];
end

%% Create clin_data_1
T=readtable(EXCEL_DIR);
clin_data_1=zeros((nr_subj),1);
for j = 1:height(T)
    for i = 1:nr_subj
        if isequal(T{j,1},SUBJECTCODES(i))
            clin_data_1(i,1)=T{j,clin_data_1_name};
            disp([SUBJECTCODES(i) ' : ' T{j,clin_data_1_name}])
        end
    end
end

%% Run multiregression
LTNP_vbm_multiregres(scanfiles, maskfiles, clin_data_1, clin_data_1_name, MASK_NAME, contrast_name,contrast_inv_name, contrast, contrast_inv,OUTPUT_DIR)

% % Define matlabbatch for making a mask with imcalc
% matlabbatch={};
% matlabbatch{1}.spm.util.imcalc.input = filelist_mask;
% matlabbatch{1}.spm.util.imcalc.output = MASK_NAME;
% matlabbatch{1}.spm.util.imcalc.outdir = {OUTPUT_DIR};
% matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)>0.3';
% matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
% 
% % Define matlabbatch for "Basic Model"
% matlabbatch{2}.spm.stats.factorial_design.dir = {OUTPUT_DIR};
% matlabbatch{2}.spm.stats.factorial_design.des.mreg.scans = filelist_scans;
% matlabbatch{2}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
% matlabbatch{2}.spm.stats.factorial_design.des.mreg.incint = 1;
% matlabbatch{2}.spm.stats.factorial_design.cov.c = gds;
% matlabbatch{2}.spm.stats.factorial_design.cov.cname = 'gds_15';
% matlabbatch{2}.spm.stats.factorial_design.cov.iCFI = 1;
% matlabbatch{2}.spm.stats.factorial_design.cov.iCC = 1;
% matlabbatch{2}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
% matlabbatch{2}.spm.stats.factorial_design.masking.tm.tm_none = 1;
% matlabbatch{2}.spm.stats.factorial_design.masking.im = 1;
% matlabbatch{2}.spm.stats.factorial_design.masking.em = {MASK_PATH};
% matlabbatch{2}.spm.stats.factorial_design.globalc.g_omit = 1;
% matlabbatch{2}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
% matlabbatch{2}.spm.stats.factorial_design.globalm.glonorm = 1;
% 
% % Define matlabbatch for "Estimate"
% matlabbatch{3}.spm.stats.fmri_est.spmmat = {[OUTPUT_DIR '/SPM.mat']};
% matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
% matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
% 
% % Define matlabbatch for "Results", then run
% matlabbatch{4}.spm.stats.con.spmmat = {[OUTPUT_DIR '/SPM.mat']};
% matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'Positive correlations';
% matlabbatch{4}.spm.stats.con.consess{1}.tcon.convec = [0 1];
% matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Negative correlations';
% matlabbatch{4}.spm.stats.con.consess{2}.tcon.convec = [0 -1];
% matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
% matlabbatch{4}.spm.stats.con.delete = 1;
% 
% % Run the job
% spm_jobman('run',matlabbatch);
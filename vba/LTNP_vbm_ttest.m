function LTNP_vbm_ttest(files_1, files_2, MASK_NAME, clin_data_1, clin_data_name_1, clin_data_2, clin_data_name_2, contrast_name,contrast_inv_name, contrast, contrast_inv, OUTPUT_DIR)

%vbm_t1_ttest(t1_files, MASK_NAME, spet_files_1, spet_files_2, clin_data_1, clin_data_name_1, clin_data_2, clin_data_name_2, OUTPUT_DIR,contrast_name,contrast_name_inv, convec, convec_inv)


% The goal of this script is to automate the VBA protocol from Leen analysis:
%       Info_VOI_Voxel_analyses.docx
%
% It will first make a GM binary mask of the average of all subjects GM
% Then make a ttest between two samples of s8wrSUVR images
% Add age as listed in a excel file as covariate
%
% Assumptions
%     - s8wrSUVR images are available
%     - wc1T1 images are also available
%     - excel file with a subjectcode per row in the first column, and the
%     corresponding age in the secund column

% Calling SPM in the outputdirectory
if isfile(MASK_NAME) 
    n=0;
    [MASK_dir,MASK_name,MASK_ext]=fileparts(MASK_NAME);
    if isequal(MASK_dir,OUTPUT_DIR)
        MASK_PATH=MASK_NAME;
    else
        MASK=[MASK_name MASK_ext];
        MASK_PATH=fullfile(OUTPUT_DIR,MASK);
        copyfile(MASK_NAME,MASK_PATH);
    end
else
    n=1;
    MASK_PATH=fullfile(OUTPUT_DIR,MASK_NAME);
    
    % Merge lists; considering these are smoothed T1 files
    sT1_files_all=[files_1;files_2];

    % Define matlabbatch for making a mask with imcalc
    matlabbatch={};
    matlabbatch{n}.spm.util.imcalc.input = sT1_files_all;
    matlabbatch{n}.spm.util.imcalc.output = MASK_NAME;
    matlabbatch{n}.spm.util.imcalc.outdir = {OUTPUT_DIR};
    matlabbatch{n}.spm.util.imcalc.expression = 'mean(X)>0.3';
    matlabbatch{n}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{n}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{n}.spm.util.imcalc.options.mask = 0;
    matlabbatch{n}.spm.util.imcalc.options.interp = 1;
    matlabbatch{n}.spm.util.imcalc.options.dtype = 4;
end

% Define matlabbatch for "Basic Model"
matlabbatch{n+1}.spm.stats.factorial_design.dir = {OUTPUT_DIR};
matlabbatch{n+1}.spm.stats.factorial_design.des.t2.scans1 = files_1;
matlabbatch{n+1}.spm.stats.factorial_design.des.t2.scans2 = files_2;
matlabbatch{n+1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{n+1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{n+1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{n+1}.spm.stats.factorial_design.des.t2.ancova = 0;
if ~isempty(clin_data_1)
    matlabbatch{n+1}.spm.stats.factorial_design.cov(1).c = clin_data_1;
    matlabbatch{n+1}.spm.stats.factorial_design.cov(1).cname = clin_data_name_1;
    matlabbatch{n+1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{n+1}.spm.stats.factorial_design.cov(1).iCC = 1;
    if ~isempty(clin_data_2)
        matlabbatch{n+1}.spm.stats.factorial_design.cov(2).c = clin_data_2;
        matlabbatch{n+1}.spm.stats.factorial_design.cov(2).cname = clin_data_name_2;
        matlabbatch{n+1}.spm.stats.factorial_design.cov(2).iCFI = 1;
        matlabbatch{n+1}.spm.stats.factorial_design.cov(2).iCC = 1;
    end
else
    if ~isempty(clin_data_2)
        matlabbatch{n+1}.spm.stats.factorial_design.cov(1).c = clin_data_2;
        matlabbatch{n+1}.spm.stats.factorial_design.cov(1).cname = clin_data_name_2;
        matlabbatch{n+1}.spm.stats.factorial_design.cov(1).iCFI = 1;
        matlabbatch{n+1}.spm.stats.factorial_design.cov(1).iCC = 1;
    end
end
matlabbatch{n+1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{n+1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{n+1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{n+1}.spm.stats.factorial_design.masking.em = {MASK_PATH};
matlabbatch{n+1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{n+1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{n+1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Define matlabbatch for "Estimate"
matlabbatch{n+2}.spm.stats.fmri_est.spmmat = {[OUTPUT_DIR '/SPM.mat']};
matlabbatch{n+2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{n+2}.spm.stats.fmri_est.method.Classical = 1;

% Define matlabbatch for "Results"
matlabbatch{n+3}.spm.stats.con.spmmat = {[OUTPUT_DIR '/SPM.mat']};
matlabbatch{n+3}.spm.stats.con.consess{1}.tcon.name = contrast_name; 
matlabbatch{n+3}.spm.stats.con.consess{1}.tcon.convec = contrast;
matlabbatch{n+3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{n+3}.spm.stats.con.consess{2}.tcon.name = contrast_inv_name;
matlabbatch{n+3}.spm.stats.con.consess{2}.tcon.convec = contrast_inv; 
matlabbatch{n+3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{n+3}.spm.stats.con.delete = 1;

% Run the job
cd(OUTPUT_DIR)
spm pet
spm_jobman('run',matlabbatch);

end
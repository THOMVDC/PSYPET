function LTNP_vbm_multiregres(files_1, maskfiles, clin_data_1, clin_data_1_name, MASK_NAME, contrast_name,contrast_inv_name, contrast, contrast_inv,OUTPUT_DIR)

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

    % Define matlabbatch for making a mask with imcalc
    matlabbatch={};
    matlabbatch{n}.spm.util.imcalc.input = maskfiles;
    matlabbatch{n}.spm.util.imcalc.output = MASK_NAME;
    matlabbatch{n}.spm.util.imcalc.outdir = {OUTPUT_DIR};
    matlabbatch{n}.spm.util.imcalc.expression = 'mean(X)>0.3';
    matlabbatch{n}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{n}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{n}.spm.util.imcalc.options.mask = 0;
    matlabbatch{n}.spm.util.imcalc.options.interp = 1;
    matlabbatch{n}.spm.util.imcalc.options.dtype = 4;
end

% Multiple regression 
matlabbatch{n+1}.spm.stats.factorial_design.dir = {OUTPUT_DIR};
matlabbatch{n+1}.spm.stats.factorial_design.des.mreg.scans = files_1;
matlabbatch{n+1}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
matlabbatch{n+1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{n+1}.spm.stats.factorial_design.cov.c = clin_data_1; % GDS
matlabbatch{n+1}.spm.stats.factorial_design.cov.cname = clin_data_1_name;
matlabbatch{n+1}.spm.stats.factorial_design.cov.iCFI = 1;
matlabbatch{n+1}.spm.stats.factorial_design.cov.iCC = 1;
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

% Define matlabbatch for "Results", then run
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



% 	Ingeven van thresholds
% 	P=0,001, Kext = lager vb 50


end

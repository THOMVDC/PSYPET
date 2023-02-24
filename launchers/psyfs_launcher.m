function psyfs(T1,PET,tracer,ref_VOI_fs,atlas_fs,atlas_fs_csv,WML,output_folder,x)

    % This version is ready to go ! -> ambition to make a single
    % dependency on spm/cat12 (not ANTS or others)
    psypet_version='v3.0'; 
    fprintf('PSYPET %s running \n', psypet_version); 
    fprintf('author: Thomas Vande Casteele \n');
    fprintf('coauthors: Patrick Dupont, Nathalie Mertens, Michel Koole \n');

%% A. BACKGROUND %%
%%%%%%%%%%%%%%%%
%
% atlas_mni and ref_VOI_mni should have same bounding box and voxelsize
% cortical thickness is estimated in the segmentation of CAT12 but only
% with freesurfer atlases (no volume based atlas), have to refer to ANTS to
% obtain a cortical thickness estimation. If ANTs is installed, uncomment
% the necessary lines below (all contain the keyword "thick")
%
% % VOIdetails should'nt contain a 0 (background), we can handle it for cat12 atlasses but
% best to avoid (cfr LTNP_VOI_stats_v4)

% Default voxelsize in LTNP_spm12_warp = 1.5; default in 

% dependencies: 
%       SPM and cat12 are in your matlab path
%       ants, afni, dcm2niix are installed and callable from terminal
%       path to scriptdir should be inserted in line 64 -> should use
%       dbstack or functioncallinfo instead

% assumptions:
%	    T1 is the path to the folder containing the T1 dcm files 
%                           OR
%       T1 is the complete path to the T1 nifti file that has been
%       processed by CAT12_7 in the same folder !
%       PET is the path to the folder containing the T1 dcm files
%       ref_VOI is the path to your ref_VOI (.nii format)
%       tracer is the name of your tracer ('UCBJ', or 'MK62', or 'FLUT')
%       output_folder is the path to your output_folder
%       !!! we assume [~,subj,~] = fileparts(output_folder) !!!
%

% atlas_mni should be a path to own atlas in nifti format, or the name of
% one of the default cat12 atlasses ('cobra','neuromorphometrics',...)
% possibility of course to add own atlas to cat12_defaults
% also obligated to put a CSV file in cat12/templates_volumes. The csv-file
% should have an header line containing the number of the ROI "ROIid", the
% abbreviation of the ROI "ROIabbr" and the full name of the ROI "ROIname".
% The GM, WM and CSF values will be extracted for all regions.

% Pet_path and seg_path should be coregistred with ants on before hand

%LTNP_Ants_coreg()

% Relies on the script of Nathalie Mertens
% Dependencies:
% To check which python installation matlab calls:
%   pyenv
% To change the python installation matlab calls (example):
%   pyenv('Version','/Library/Frameworks/Python.framework/Versions/3.6/bin/python3.6')

% On iPsychiater, the directory where python rbv script is located is 
% script_dir='/Volumes/LaCie/Thomas/Projects/L3D/SCRIPTS/NATHALIE_MERTENS/'
% script_name='PVC_RBV_calculation_Nathalie.py'


% TO DO:
%
%   1/ center, crop (with LTNP_crop.m, the best one) , acpc T1  (matlab function, also callable by python for
%   FLAIR, CUBE, etc)
%   2/ realign and mean PET, then center, crop and acpc
%   3/ coregister PET 2 T1
%   4/ segment T1
%   5/ create suvr in mni space
%   6/ bring to patient space
%   7/ write a SUVR (patient space), rSUVR (coregistred to the mr) and
%   wSUVR (warped to mni)
%
% !!!!!!!!!  antsreg seems to work better than spmreg for FLUT data atleast  !!!!
% !!!!!!!!!  T1 and PET should be in nifti (dcm 2 suv !) or processed by
% python with dcm 2 suv ?????? See function "LTNP_dcm2nii_SUV"
% 

%% B. SETTINGS %%
%%%%%%%%%%%%%%

if nargin<9
    use_defaults=1;
else
    use_defaults=0;
end

if use_defaults==1
    
    % Grab CAT version software
    %CAT_version=cat_version;
    CAT_version='fastsurfer'; 
    
    % Defining processing_id
    script = mfilename('fullpath');
    if isequal(script,[])
        error('Script not found')
    else
        fprintf('PSYPET script located in %s \n', script);         
    end
    script_ext='.m';
    [script_folder,script_name,~]=fileparts(script);
    rbv_script_dir=script_folder;
    scriptid= [script_name '_' psypet_version];
    processing_id = char([scriptid '_processed_' CAT_version '_' date]);
    
    % Defining directories
    [main_dir,subjectcode,~]    = fileparts(output_folder); % directory where the folders of each subject can  be found  
            
    % Define tracer and tracer specific info
    [~,ref_VOI_name,~] = fileparts(ref_VOI_mni);
    
    % Define current subject 
    subjdir = fullfile(output_folder, processing_id);
    cd(output_folder)
    mkdir(processing_id)
    
    % Define PET specific directories
    subjdir_PET = fullfile(subjdir,tracer); 
    subjdir_PETdcm = fullfile(subjdir_PET,'DCM');
    subjdir_PETnii = fullfile(subjdir_PET,'NIFTI');
    %subjdir_PETpvcmg = fullfile(subjdir_PET,'PVC_MG');
    subjdir_PETpvcrbv = fullfile(subjdir_PET,'PVC_RBV');
    %subjdir_PETw = fullfile(subjdir_PET,'WARPED');
    %subjdir_PETcm = fullfile(subjdir_PET,'COLORMAPS');
    subjdir_PETcor= fullfile(subjdir_PET,'COREG');
    subjdir_PETsuvr= fullfile(subjdir_PET,'SUVR');
    subjdir_PETreal= fullfile(subjdir_PET,'REALIGNED');
    mkdir(subjdir_PETreal)
    mkdir(subjdir_PETsuvr)
    %mkdir(subjdir_PETcm)
    mkdir(subjdir_PETdcm)
    mkdir(subjdir_PETnii)
    %mkdir(subjdir_PETpvcmg)
    mkdir(subjdir_PETpvcrbv)
    %mkdir(subjdir_PETw)
    mkdir(subjdir_PETcor)
    copyfile(PET, subjdir_PETdcm)
    
    % Define T1 specific directories
    if endsWith(T1,'.nii')
        [subjdir_ANAT,~,~]=fileparts(T1);
        subjdir_ANATsegm=fullfile(subjdir_ANAT,'CAT12');
        subjdir_ANATatlas=fullfile(subjdir_ANAT,'RBV_ATLASSES');
        subjdir_ANATmask=fullfile(subjdir_ANAT,'MASKS');
        %subjdir_ANATthick=fullfile(subjdir_ANAT,'THICKNESS');
    else
        subjdir_ANAT = fullfile(subjdir,'ANAT'); 
        subjdir_ANATsegm = fullfile(subjdir_ANAT,'CAT12');
        subjdir_ANATatlas= fullfile(subjdir_ANAT,'RBV_ATLASSES');
        subjdir_ANATmask=fullfile(subjdir_ANAT,'MASKS');
        %subjdir_ANATthick=fullfile(subjdir_ANAT,'THICKNESS');
        mkdir(subjdir_ANATsegm)
        %mkdir(subjdir_ANATthick)
        mkdir(subjdir_ANATatlas)
        mkdir(subjdir_ANATmask)
    end

        
    % Define tracer specific values
    if strcmp(tracer,'FLUT')
        nr_frames=6;
        SUVR_start_time    = 90; % in min
        SUVR_end_time      = 120; % in min
    elseif strcmp(tracer,'UCBJ')
        nr_frames=6; %8
        SUVR_start_time    = 60; % in min
        SUVR_end_time      = 90; % in min
    elseif strcmp(tracer,'MK62')
        nr_frames=6;
        SUVR_start_time    = 90; % in min
        SUVR_end_time      = 120; % in min
    elseif strcmp(tracer,'PIB')
        nr_frames=4;
        SUVR_start_time    = 90; % in min
        SUVR_end_time      = 120; % in min
    else
        nr_frames='unknown';
        SUVR_start_time    = []; % in min
        SUVR_end_time      = []; % in min
    end
    
    % Defining SUVR parameters
%     additional_smooth_parametric = 0; % kernel size in mm; isotropic Gaussian 3D smoothing
%     GM_threshold                 = 0.3;
%     WM_threshold                 = 0.2;
    threshold_totmove = 1; % maximum scan to scan movement in mm (used in function LCN12_analyze_headmovement)
    threshold_translation = 3; % maximum overall movement in mm (used in function LCN12_analyze_headmovement)
    threshold_rotation    = 3; % maximum overall rotation in degrees (used in function LCN12_analyze_headmovement)

    % Define PVC parameters
    %FWHM = 5; % FHWM of the image resolution (assumed isotropic) in mm
    
    % Define if you want to output white matter hyperintensities
%    WMHC=1;

    % Atlasses you can choose from inside cat12
%     ATLASSES_cat12={ % in patient space, those are generated during cat12 segmentation, you can add any you can find in path_to_VOIdetails
%         'hammers'
%         'cobra'
%         'neuromorphometrics'
%     };
%     cat_dir = which('cat12');
%     cat_dir = cat_dir(1:end-8);
%     %path_to_VOIdetails=fullfile(cat_dir,'templates_MNI152NLin2009cAsym');
%     path_to_VOIdetails=fullfile(cat_dir,'templates_volumes');
    
    %
    %
    % Grab templates
    %mni_origin=fullfile(script_folder,'/templates/ACsphere_bin.nii.gz');
    %mni=fullfile(script_folder,'/templates/Template_T1_IXI555_MNI152_GS.nii');

       
elseif use_defaults==0

%% Use of configfile

%      % Reading json config file
      config=jsondecode(fileread(x));
% % 
%      % Defining directories
%      main_dir    = config.directories.output_dir; % directory where the folders of each subject can  be found
%      script_dir  = config.directories.L3D_CORE_script_dir;
% %     
% %     % Defining subjects and number of subjects
%       SUBJECTSCODES=config.subjects;
%       SUBJECTSDIRS=cell(size(SUBJECTSCODES));
%       for subj=1:length(SUBJECTSCODES)
%           SUBJECTSDIRS(subj)=fullfile(main_dir,SUBJECTSCODES(subj));
%       end
       subjectcode= config.subject;
%     subjectdir = char(SUBJECTSDIRS(subj,:));  
%     
%     % Defining tracer and tracer specific info
%     if config.analysis.UCBJ.execute=='yes'
%         tracer='UCBJ';
%         ref_VOI = fullfile(script_dir,'Templates','Mask_VOI_SO_AtlasspacePMOD.nii');
%         ref_VOI_name = 'semi_ovale';
%         %infostring_meta = 'UCBJ_frames';
%         SUVR_start_time    = 60; % in min
%         SUVR_end_time      = 90; % in min
%         if strcmp(config.analysis.UCBJ.LTNP.PVC,'yes')
%             PVC='PVC';
%         else
%             PVC='';
%         end
%     elseif config.analysis.MK62.execute=='yes'
%         tracer='MK62';
%         ref_VOI = fullfile(script_dir,'Templates','Mask_VOI_iCBL_AtlasspacePMOD.nii');
%         ref_VOI_name = 'iCBL';
%         %infostring_meta = 'MK62_frames';
%         SUVR_start_time    = 90; % in min
%         SUVR_end_time      = 120; % in min
%         if strcmp(config.analysis.MK62.LTNP.PVC,'yes')
%             PVC='PVC';
%         else
%             PVC='';
%         end
%     elseif config.analysis.FLUT.execute=='yes'
%         tracer='FLUT';
%         ref_VOI = fullfile(script_dir,'Templates','resl_to_wmc1_SPM12_aal_whole_cerebellum.img');
%         ref_VOI_name = 'cerebellum';
%         %infostring_meta = 'FLUT_frames';
%         SUVR_start_time    = 90; % in min
%         SUVR_end_time      = 120; % in min
%         if strcmp(config.analysis.FLUT.LTNP.PVC,'yes')
%             PVC='PVC';
%         else
%             PVC='';
%         end
%     else
%         fprintf('This script is not implemented for the specified tracer');
%         return
%     end
%     
%     % Defining processing_id
%     functioncallinfo = dbstack;
%     scriptname = functioncallinfo.name;
%     processing_id = char([scriptname '_processed_' config.directories.processing_ID]);
%     
% %     % Define subjectspecific directories
%       subjdir_ANATsegm = fullfile(config.directories.datadir_SEGM,subjectcode);
% %     %subjdir_frame_definition= fullfile(subjectdir,tracer,config.directories.datadir_frame_definition);
% %     subjdir_PET = fullfile(subjectdir,tracer,config.directories.datadir_PET); 
%       subjdir_ANAT = fullfile(subjectdir,'T1',config.directories.datadir_ANAT); 
%     subjdir_PETprocessed = fullfile(subjectdir,tracer,processing_id);

%     % Defining masks and VOIs
%     BRAINmask_path = 
%     GM_mask_file = 
% 
%     % Defining SUVR parameters
%     additional_smooth_parametric = 
%     GM_threshold                 = 
%     threshold_totmove =
%     threshold_translation = 
%     threshold_rotation    = 
%     % data_space                   = 
% 
%     % Define PVC parameters
%     FWHM = 
% 
%     % Choose atlasses, calculate number of atlasses
%     path_to_atlasses=
%     ATLASSES={
%         'Hammers_N30R83_1mm'
%     };
%     nr_atlas = length(ATLASSES);

end

%% C. PROCESSING %%
%%%%%%%%%%%%%%%% 

%%      1/ Copy settings variables and create logfile
% ---------------------------------

% Copy the current script in subjdir_PETprocessed
copyfile(strcat(script,script_ext),fullfile(subjdir,[script_name '.txt']))

% Save the current variables in subjdir_PETprocessed 
cd(subjdir)
save(['settings_variables_' tracer '_' subjectcode '_' processing_id '.mat'])
list_settings_variables=who;

% Create a logfile
fprintf('working on subject %s for tracer %s \n',subjectcode,tracer);   
name_logfile = fullfile(subjdir,char([processing_id '_' subjectcode '.txt'])); 
fid  = fopen(name_logfile,'a+'); 

% Talk to it, tell the settings that are used
fprintf(fid,'subject = %s\n',subjdir);
fprintf(fid, '%s\n',script_name);
fprintf(fid,'%c','-'*ones(1,30));
fprintf(fid,'\n');
tmp  = clock;  
fprintf(fid,'Processing started: %s at %i h %i min %i s\n',date,tmp(4), tmp(5),round(tmp(6)));
fprintf(fid,'\n');
fprintf(fid,'Settings\n');
fprintf(fid,'reference_VOI = %s\n',ref_VOI_mni);
fprintf(fid,'additional_smooth_parametric = %i mm \n',additional_smooth_parametric);
fprintf(fid,'GM_threshold = %2.1f\n',GM_threshold);
fprintf(fid,'threshold (scan to scan movement in mm) = %4.2f \n',threshold_totmove);
fprintf(fid,'threshold (overall translation in mm)   = %4.2f \n',threshold_translation);
fprintf(fid,'threshold (overall rotation in degrees) = %4.2f \n',threshold_rotation);

%%      2/ Preprocess T1 (dcm2nii, crop, center)
% --------------------------------------

% Grab atlasses
% if any(strcmp(ATLASSES_cat12,atlas_mni))
%     atlas_mni_name=atlas_mni;
%     VOIdetails_path=fullfile(path_to_VOIdetails,[atlas_mni '.csv']);
% else
%     [~, atlas_mni_name, ~]=fileparts(atlas_mni);
%     VOIdetails_path=atlas_csv;
% end

% Grab voxelsize
voxelsize=LTNP_get_voxelsize(ref_VOI_fs);

% Check if T1 needs to be processed
if endsWith(T1,'.nii')
    
    % T1 doesn't need to be processed
    GoT1=false;
    
    % Talk to log
    fprintf(fid,'T1 preprocessing already done \n');
    accT1nii=T1;
    
    % We assume it's the right subject (cfr coregistration)
    ptname_T1='unknown';
    ptID_T1='unknown';
    ptname_PET='unknown';
    ptID_PET='unknown';
    
    % Grab ref image name and ext
    [~,ref_image_name,ref_image_ext]=fileparts(accT1nii);
    
    % Set paths we'll need
%     deformation_field = fullfile(subjdir_ANATsegm,'mri',['y_' ref_image_name ref_image_ext]);
%     invdeformation_field = fullfile(subjdir_ANATsegm,'mri',['iy_' ref_image_name ref_image_ext]);
%     CAT12_vol_thick=load(fullfile(subjdir_ANATsegm,'report',['cat_' ref_image_name '.mat']));
%     str=CAT12_vol_thick.S.catlog{end-4, 1};
%     IQR=extractBetween(str,'(IQR):  ','%');
%     cat12_mean_thickness=cell2mat(CAT12_vol_thick.S.subjectmeasures.dist_thickness);
%     TIV=CAT12_vol_thick.S.subjectmeasures.vol_TIV;
%     %thickness_path=fullfile(subjdir_ANATthick,['thickness_p1' ref_image_name ref_image_ext]); 
%     %wthickness_path=fullfile(subjdir_ANATthick,['thickness_mwp1' ref_image_name ref_image_ext]);    
%     rbv_atlas_path=fullfile(subjdir_ANATatlas,['rbv_segm_' atlas_mni_name '_' ref_image_name ref_image_ext]);   
    BRAINmask_path=fullfile(subjdir_ANATmask,['GMWMCSF_mask_p' ref_image_name ref_image_ext]);
    wBRAINmask_path=fullfile(subjdir_ANATmask,['GMWMCSF_mask_mwp' ref_image_name ref_image_ext]);
    GMmask_path=fullfile(subjdir_ANATmask,['mask_p1' ref_image_name ref_image_ext]);
    WMmask_path=fullfile(subjdir_ANATmask,['mask_p2' ref_image_name ref_image_ext]);
    CSFmask_path=fullfile(subjdir_ANATmask,['mask_p3' ref_image_name ref_image_ext]);
    WMHmask_path=fullfile(subjdir_ANATmask,['mask_p7' ref_image_name ref_image_ext]);
%     wGMmask_path=fullfile(subjdir_ANATmask,['mask_mwp1' ref_image_name ref_image_ext]);
%     wWMmask_path=fullfile(subjdir_ANATmask,['mask_mwp2' ref_image_name ref_image_ext]);
%     wCSFmask_path=fullfile(subjdir_ANATmask,['mask_mwp3' ref_image_name ref_image_ext]);
%     wWMHmask_path=fullfile(subjdir_ANATmask,['mask_mwp7' ref_image_name ref_image_ext]);

else
    
    % T1 has to be processed
    GoT1=true;
    
    % Talk to log
    fprintf(fid,'T1 processing started');
    
    % Run
    T1nii_name=['T1_' subjectcode '.nii'];
    [accT1nii, ptname_T1, ptID_T1]=LTNP_preproc_T1_spm(T1,subjdir_ANAT,T1nii_name);
   
    % Talk to log
    fprintf(fid,'T1 Participant: %s\n',ptname_T1);
    fprintf(fid,'T1 Participant ID: %s\n',ptID_T1);
    fprintf(fid,'T1 processing finished \n');
  
end

%%      3/ Segment T1
%       -------------

if GoT1
    
    % Grab name
    [~,ref_image_name,ref_image_ext]=fileparts(accT1nii);
    
    % Talk to logfile
    fprintf(fid,'T1 segmentation started');

    % Segment
    [IQR, TIV, cat12_mean_thickness, deformation_field, invdeformation_field]=LTNP_cat12_7_segment(accT1nii, subjdir_ANATsegm, WMHC, atlas_mni, voxelsize); % also creates the deformation field (warping parameters) to MNI

    % Talk to log file
    tmp  = clock;
    fprintf(fid,'Segmentation ended at: %i h %i min %i s\n',tmp(4), tmp(5),round(tmp(6)));
    fprintf(fid,'Image Quality Rating: s% s\n',char(IQR));
    fprintf(fid,'T1 segmentation calculation finished \n');
    
end

% Grab output maps generated from segmentation (you'll need it later)
%GWC_path=fullfile(subjdir_ANATsegm,'mri',['p0' ref_image_name ref_image_ext]); % GM+WM+CSF in patientspace
GM_path=fullfile(subjdir_ANATsegm,'mri',['p1' ref_image_name ref_image_ext]); % GM in patientspace
WM_path=fullfile(subjdir_ANATsegm,'mri',['p2' ref_image_name ref_image_ext]); % WM in patientspace
CSF_path=fullfile(subjdir_ANATsegm,'mri',['p3' ref_image_name ref_image_ext]); % CSF in patientspace
WMH_path=fullfile(subjdir_ANATsegm,'mri',['p7' ref_image_name ref_image_ext]); % White matter hyperintensities in patientspace
%wGWC_path=fullfile(subjdir_ANATsegm,'mri',['mwp0' ref_image_name ref_image_ext]); % GM+WM+CSF in patientspace
wGM_path=fullfile(subjdir_ANATsegm,'mri',['mwp1' ref_image_name ref_image_ext]); % GM in mni
wWM_path=fullfile(subjdir_ANATsegm,'mri',['mwp2' ref_image_name ref_image_ext]); % WM in mni
wCSF_path=fullfile(subjdir_ANATsegm,'mri',['mwp3' ref_image_name ref_image_ext]); % CSF in mni
wWMH_path=fullfile(subjdir_ANATsegm,'mri',['mwp7' ref_image_name ref_image_ext]); % WMH in mni

% Grab output atlas generated from segmentation (you'll need it later)
path_to_atlasses_cat12=fullfile(subjdir_ANATsegm,'mri_atlas'); % Grab atlasses in patient space
atlas_path=fullfile(path_to_atlasses_cat12,[atlas_mni_name '_' ref_image_name ref_image_ext]);
watlas_path=fullfile(path_to_atlasses_cat12,['w' atlas_mni_name '_' ref_image_name ref_image_ext]);

% Grab voxelsize of MR patient space
voxelsize_pt=LTNP_get_voxelsize(accT1nii);

%%      4/ Create masks
% ---------------------------

if GoT1

    % GM, WM and CSF masks
    [GMmask_path,WMmask_path,CSFmask_path,BRAINmask_path,~]=LTNP_make_labelimage(GM_path,WM_path,CSF_path,subjdir_ANATmask); % Creates the masks
    
    % wGM, wWM and wCSF masks
    [wGMmask_path,wWMmask_path,wCSFmask_path,wBRAINmask_path,~]=LTNP_make_labelimage(wGM_path,wWM_path,wCSF_path,subjdir_ANATmask);

    % WMH masks
    WMHmask_path=fullfile(subjdir_ANATmask,['mask_p7' ref_image_name ref_image_ext]);
    wWMHmask_path=fullfile(subjdir_ANATmask,['mask_mwp7' ref_image_name ref_image_ext]);
    if WMHC==1
        LTNP_binarize(WMH_path,WMHmask_path, 0, Inf);
        LTNP_binarize(wWMH_path,wWMHmask_path, 0, Inf);
    end
end

%%      5/ Generate a cortical thickness map
% if GoT1
%     
%     % Talk to logfile
%     fprintf('subject %s: Generating thickness map for subject: \n',subjectcode);
%     
%     % The use of this is largely justified by A large-scale comparison of cortical thickness and volume methods for measuring Alzheimer's disease severity
%     thickness_path=LTNP_ANTs_thickness(GM_path,WM_path,CSF_path,subjdir_ANATthick);
%     wthickness_path=LTNP_ANTs_thickness(wGM_path,wWM_path,wCSF_path,subjdir_ANATthick);
% end

%%      6/ MAKE a rbv atlas map

% Create rbv segmentation and brain mask
if GoT1
    [rbv_atlas_path,~,~]=LTNP_create_rbv_segm_v13(GMmask_path, WMmask_path, CSFmask_path, atlas_path, subjdir_ANATatlas);
    %[wrbv_atlas_path,~,~]=LTNP_create_rbv_segm_v13(wGM_path, wWM_path, wCSF_path, watlas_path, subjdir_ANATatlas);
end

%%      7/ Make T1 statistics

% Make segmentation stats patient space
[GM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,GM_path,GMmask_path);
[WM_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,WM_path,WMmask_path);
[CSF_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,CSF_path,CSFmask_path);
if WMHC==1
    [WMH_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,WMHmask_path);
else
    WMH_table=[];
end
stats=horzcat(GM_table,WM_table(:,2:end),CSF_table(:,2:end),WMH_table(:,2:end));

% Make segmentation stats mni space
[wGM_table,~,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,wGM_path);
[wWM_table,~,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,wWM_path);
[wCSF_table,~,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,wCSF_path);
if WMHC==1
    [wWMH_table,~,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,wWMH_path);
else
    wWMH_table=[];
end
wstats=horzcat(wGM_table,wWM_table(:,2:end),wCSF_table(:,2:end),wWMH_table(:,2:end));

% Get stats out of thickness
%[thickness_table,~,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,thickness_path);
%[wthickness_table,~,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,wthickness_path);

% Concatenate with segmentation stats
%stats=horzcat(stats,thickness_table(:,2:end));
%wstats=horzcat(wstats,wthickness_table(:,2:end));

%%      8/ Convert PET to nii, and transform to SUV frames
% --------------------------------------------

if endsWith(PET,'.nii')   
    
    % Talk to log
    fprintf(fid,'We assume PET preprocessing already done \n');
    
    % We assume it's the right person
    ptname_PET=[ptname_PET ''];
        
else
    
    % Talk to log
    fprintf(fid,'Centering and SUV calculation started');
       
    % Run SUV calculation
    SUVname=['SUV_' tracer '_' subjectcode];
    [ptname_PET,ptID_PET, bodyweight, injdosis, acqtime, injtime, halftime]=LTNP_dcm2SUV(subjdir_PETdcm,subjdir_PETnii,SUVname);   
    
    % Talk to log
    fprintf(fid,'PET Participant: %i g s\n',ptname_PET);
    fprintf(fid,'PET Participant ID: %i g s\n',ptID_PET);
    fprintf(fid,'Bodyweight: %i g s\n',bodyweight);
    fprintf(fid,'Injection dosis: %i bq s\n',injdosis);
    fprintf(fid,'Acquisition time: %i s\n',acqtime);
    fprintf(fid,'Injection time: %i s\n',injtime);
    fprintf(fid,'Halftime: %i s s\n',halftime);
    fprintf(fid,'SUV calculation finished \n');
    
end

%%      9/ Realign PET frames and summate
% -------------------------------

if endsWith(PET,'.nii') 
    
    % Talk to log
    fprintf(fid,'Realignment of PET frames already done n\');
    
else
    
    % Talk to log
    fprintf(fid,'Realignment of PET frames started');
    
    % Run if nr_frames corresponds to the number of nifti files found
    filelist_SUV= dir(fullfile(subjdir_PETnii,'*.nii'));
    if nr_frames==length(filelist_SUV)
        [rSUVnii_path]=LTNP_spm12_realign_and_summate(filelist_SUV,subjdir_PETreal);
    else
        error('Number of frames found in %i does not correspond to the number of frames expected %i s\n',length(filelist_SUV),nr_frames)
    end
    
    % Talk to log
    fprintf(fid,'Realignment finished \n');
    
end

%%      10/ Coregister PET to T1
%       ----------------------

% Check if participant PET is corresponding to participant T1
ptname_T1(isspace(ptname_T1)) = []; % remove spaces
ptname_T1=lower(ptname_T1); % lower all letters
ptname_T1=replace(ptname_T1,'_',''); % remove underscores
ptname_PET(isspace(ptname_PET)) = []; % remove spaces
ptname_PET=lower(ptname_PET); % lower all letters
ptname_PET=replace(ptname_PET,'_',''); % remove underscores

if strcmp(ptname_T1,ptname_PET) || strcmp(ptID_T1,ptID_PET) || strcmp(ptname_T1,ptID_PET) || strcmp(ptname_PET,ptID_T1)
   fprintf(fid,'T1 and PET corresponding');
else
   fprintf(fid,'T1 and PET not corresponding');
end

% Talk to config file
fprintf(fid,'\n');
fprintf(fid,'Coregistration started \n');

% First center PET image (T1 already aligned on the AC)
LTNP_center(rSUVnii_path)

% Let the magic happen
%rigid='r';
rrrSUVnii=LTNP_spm12_coregister_reslice(accT1nii, rSUVnii_path, subjdir_PETcor); % first rigid coreg with spm
%rrrSUVnii=LTNP_ANTs_coregister(accT1nii,rrSUVnii,rigid, subjdir_PETcor);  % now second coreg with ants

% Talk to config file
fprintf(fid,'Coregistration ended \n');

%%      11/ Apply MG PVC on PET
%        --------------------

% Talk to log
fprintf(fid,'Applying Muller Gartner partial volume correction');

% Execute PVC_MG
[rrrSUVnii_pvc_MGorig_path,rrrSUVnii_pvc_MGmodi_path,mean_WM_value] = LTNP_PVC_MG(rrrSUVnii,GM_path,WM_path,FWHM,subjdir_PETpvcmg);

% Talk to log file
fprintf(fid,'\n');
fprintf(fid,'the mean WM value for correction is %i \n',mean_WM_value);
fprintf(fid,'PVC finished \n');

%%      12/ Apply RBV PVC on PET
        % ---------------------

% Talk to log
fprintf(fid,'Applying Region Based partial volume correction');
        
% Execute RBV
cd(subjdir_PETpvcrbv) % change directory to output directory for RBV script
[rrrSUVnii_pvc_RBV_path]=LTNP_PVC_RBV(subjectcode,rbv_script_dir,rrrSUVnii,rbv_atlas_path,subjdir_PETpvcrbv);

% Move files to output_folder
mov=dir([subjdir_PET '/PVC*.*']);
for d=1:numel(mov)
    movefile(fullfile(subjdir_PET,mov(d).name),subjdir_PETpvcrbv);
end
[~,rrrSUVnii_pvc_RBV_path_img,rrrSUVnii_pvc_RBV_path_ext]=fileparts(rrrSUVnii_pvc_RBV_path);
rrrSUVnii_pvc_RBV_path=fullfile(subjdir_PETpvcrbv,[rrrSUVnii_pvc_RBV_path_img rrrSUVnii_pvc_RBV_path_ext]);

% Talk to log
fprintf(fid,'Region Based partial volume correction done');

%%      13/ Apply the warping parameters to the PET data & ref mask
%        -------------------------------------------

% Talk to log file
fprintf(fid,'\n');
tmp  = clock;
fprintf(fid,'Warping PET started: %s at %i h %i min %i s\n',date,tmp(4), tmp(5),round(tmp(6)));

% Bring PET data to mni
wrrrSUVnii=LTNP_spm12_warp(rrrSUVnii,deformation_field,subjdir_PETw,voxelsize);
wrrrSUVnii_pvc_MGorig_path=LTNP_spm12_warp(rrrSUVnii_pvc_MGorig_path,deformation_field,subjdir_PETw,voxelsize);
wrrrSUVnii_pvc_MGmodi_path=LTNP_spm12_warp(rrrSUVnii_pvc_MGmodi_path,deformation_field,subjdir_PETw,voxelsize);
wrrrSUVnii_pvc_RBV_path=LTNP_spm12_warp(rrrSUVnii_pvc_RBV_path,deformation_field,subjdir_PETw,voxelsize);

% Bring RBV atlas to mni
wrbv_atlas_path=LTNP_spm12_warp(rbv_atlas_path,deformation_field,subjdir_PETpvcrbv,voxelsize);

%Bring ref_VOI to patient space
[ref_VOI_patient]=LTNP_spm12_warp_ROI(ref_VOI_mni,invdeformation_field,subjdir_PETw,voxelsize_pt);

% Talk to log file
tmp  = clock;
fprintf(fid,'Warping PET ended at: %i h %i min %i s\n',tmp(4), tmp(5),round(tmp(6)));

% Clean warped directory from the orginal images that have been copied
del=dir([subjdir_PETw '/*.nii']);
del=del(~startsWith({del.name}, 'w'));
for d=1:numel(del)
    delete(del(d).name);
end

%%      14/ Create an output data cell with PET, MR, and atlas data
% --------------------------------------------------------------------

% Outdata
PETMR_outdata=cell(8,12);
PETMR_outdata{1,1}=rrrSUVnii;
PETMR_outdata{2,1}=rrrSUVnii_pvc_MGorig_path;
PETMR_outdata{3,1}=rrrSUVnii_pvc_MGmodi_path;
PETMR_outdata{4,1}=rrrSUVnii_pvc_RBV_path;
PETMR_outdata{5,1}=wrrrSUVnii;
PETMR_outdata{6,1}=wrrrSUVnii_pvc_MGorig_path;
PETMR_outdata{7,1}=wrrrSUVnii_pvc_MGmodi_path;
PETMR_outdata{8,1}=wrrrSUVnii_pvc_RBV_path;

PETMR_outdata{1,4}='SUVR';
PETMR_outdata{2,4}='pvc_MGorig_SUVR';
PETMR_outdata{3,4}='pvc_MGmodi_SUVR';
PETMR_outdata{4,4}='pvc_RBV_SUVR';
PETMR_outdata{5,4}='wSUVR';
PETMR_outdata{6,4}='wpvc_MGorig_SUVR';
PETMR_outdata{7,4}='wpvc_MGmodi_SUVR';
PETMR_outdata{8,4}='wpvc_RBV_SUVR';

% Fill already 5th,6th,and 7thcolumn
for item=1:8
    PETMR_outdata{item,5}=IQR;
    PETMR_outdata{item,6}=mean_WM_value;
    PETMR_outdata{item,7}=TIV;
    PETMR_outdata{item,8}=cat12_mean_thickness;
    PETMR_outdata{item,11}=bodyweight;
    PETMR_outdata{item,12}=injdosis;
end

%%      15/ Calculate SUVR and ROI data
%       --------------------------

% Talk to logfile
tmp  = clock;
fprintf(fid,'\n');
fprintf(fid,'Calculating SUVR image started: %s at %i h %i min %i s\n',date,tmp(4), tmp(5),round(tmp(6)));
fprintf('\n');
fprintf('reading data of subject %s ... \n',subjectcode);

for k=1:8
       
    % Read data
    [SUV, Vref]=LCN12_read_image(PETMR_outdata{k,1});
    if k<5 % we'll read patient space data
        GMimg=LCN12_read_image(GM_path,Vref);
        WMimg=LCN12_read_image(WM_path,Vref);
        WMH=WMH_path;
        refVOIimg=LCN12_read_image(ref_VOI_patient,Vref);
        brain_mask=LCN12_read_image(BRAINmask_path,Vref); 
    else % we'll read mni space data
        GMimg=LCN12_read_image(wGM_path,Vref);
        WMimg=LCN12_read_image(wWM_path,Vref);
        WMH=wWMH_path;
        refVOIimg=LCN12_read_image(ref_VOI_mni,Vref);
        brain_mask=LCN12_read_image(wBRAINmask_path,Vref);
        atlas_path=watlas_path;
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
        fprintf(fid, '\n');
        fprintf(fid, 'SUVR calculated for FLUT. ');
        [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_flut(SUV,refVOI_thresholded,GMimg,GM_threshold,brain_mask_thresholded);  
    elseif strcmp(tracer,'UCBJ')
        fprintf(fid, '\n');
        fprintf(fid, 'SUVR calculated for UCBJ. ');
        [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_ucbj(SUV,refVOI_thresholded,WMimg,WM_threshold,brain_mask_thresholded,WMHimg);
    elseif strcmp(tracer,'MK62')
        fprintf(fid, '\n');
        fprintf(fid, 'SUVR calculated for MK62. ');
        [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_mk62(SUV,refVOI_thresholded,GMimg,GM_threshold,brain_mask_thresholded);
    elseif strcmp(tracer,'PIB')
        fprintf(fid, '\n');
        fprintf(fid, 'SUVR calculated for Pib. ');
        [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_pib(SUV,refVOI_thresholded,GMimg,GM_threshold,brain_mask_thresholded);
    else
        fprintf(fid, '\n');
        fprintf(fid, 'No tracer specific processing defined. SUVR calculated in a standard way. ');
        [SUVR_img, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR(PET,refVOI_thresholded,brain_mask_thresholded);
    end
        
    % Save SUVR, ref mask, nr_voxels_refVOI and ref_value
    PETMR_outdata{k,2}=nr_voxels_refVOI;
    PETMR_outdata{k,3}=ref_value;
    [~, SUV_name, ~]=fileparts(PETMR_outdata{k,1});
    outputname_refmask = fullfile(subjdir_PETsuvr,['ref_mask_SUVR_' SUV_name '_' ref_VOI_name '_' tracer '_' subjectcode '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min_refVOI.nii']);
    SUVR_path = fullfile(subjdir_PETsuvr,['SUVR_' SUV_name '_' ref_VOI_name '_' tracer '_' subjectcode '_' num2str(SUVR_start_time) 'min_' num2str(SUVR_end_time) 'min.nii']); 
    LCN12_write_image(SUVR_img,SUVR_path,'SUVR',Vref.dt(1),Vref); 
    LCN12_write_image(ref_mask,outputname_refmask,'SUVR_refVOI',Vref.dt(1),Vref);
    
    % Talk to logfile
    tmp  = clock;
    fprintf(fid,'\n');
    fprintf(fid,'SUVR_start_time = %i min \n',SUVR_start_time);
    fprintf(fid,'SUVR_end_time = %i min \n',SUVR_end_time);
    fprintf(fid,'name of the reference VOI %s\n',ref_VOI_name);
    fprintf(fid,'Reference value: %i \n',ref_value);
    fprintf(fid,'number of voxels in the reference VOI = %i\n',nr_voxels_refVOI);
    fprintf(fid,'Calculating SUVR image ended at: %i h %i min %i s\n',tmp(4), tmp(5),round(tmp(6)));
    fprintf('\n');
    fprintf('Calculation of SUVR image %i from subject %s ended ... \n',k,subjectcode);
    
    %% Make PET statistics
    % ----------------
       
    % Make a csv file of the extra 3 ROIs from the rbv segmentation (cfr. LTNP_create_rbv_segm.m)
    [rbv_atlas,~]=LCN12_read_image(rbv_atlas_path);
    VOImax=unique(rbv_atlas);
    ROIid=[VOImax(end-2); VOImax(end-1); VOImax(end)];
    ROIabbr={'WM';'CSF';'MENGSKULL'};
    ROIname={'White Matter';'Cerebrospinal fluid'; 'Meninges and skull'};
    table_VOIdetails_rbv=table(ROIid,ROIabbr,ROIname);
    VOIdetails_rbv=fullfile(subjdir_ANATatlas,'VOIdetails_rbv.csv');
    writetable(table_VOIdetails_rbv,VOIdetails_rbv);
    
    % Let it run for atlas_path
    if k<5
        [SUVR_GM_table,GM_colormap,Vrefcolormap]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,GMmask_path);
        [SUVR_WM_table,WM_colormap,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,WMmask_path);
        [SUVR_CSF_table,CSF_colormap,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,CSFmask_path);
        if WMHC==1
            [SUVR_WMH_table,WMH_colormap,~]=LTNP_VOI_stats_v7(atlas_path,VOIdetails_path,SUVR_path,WMHmask_path);
        else
            SUVR_WMH_table=[];
            WMH_colormap=zeros(size(GMimg));
        end
        [SUVR_MENG_table,~,~]=LTNP_VOI_stats_v7(rbv_atlas_path,VOIdetails_rbv,SUVR_path);
        PETMR_outdata{k,9}=horzcat(stats,SUVR_GM_table(:,2:end),SUVR_WM_table(:,2:end),SUVR_CSF_table(:,2:end),SUVR_WMH_table(:,2:end));
        PETMR_outdata{k,10}=horzcat(PETMR_outdata{k,10},SUVR_MENG_table(:,2:end));
    else
        [wSUVR_GM_table,wGM_colormap,wVrefcolormap]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,SUVR_path,wGMmask_path);
        [wSUVR_WM_table,wWM_colormap,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,SUVR_path,wWMmask_path);
        [wSUVR_CSF_table,wCSF_colormap,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,SUVR_path,wCSFmask_path);
        if WMHC==1
            [wSUVR_WMH_table,wWMH_colormap,~]=LTNP_VOI_stats_v7(watlas_path,VOIdetails_path,SUVR_path,wWMHmask_path);
        else
            wSUVR_WMH_table=[];
            wWMH_colormap=zeros(size(GMimg));
        end
        [wSUVR_MENG_table,~,~]=LTNP_VOI_stats_v7(wrbv_atlas_path,VOIdetails_rbv,SUVR_path);
        PETMR_outdata{k,9}=horzcat(wstats,wSUVR_GM_table(:,2:end),wSUVR_WM_table(:,2:end),wSUVR_CSF_table(:,2:end),wSUVR_WMH_table(:,2:end));
        PETMR_outdata{k,10}=horzcat(PETMR_outdata{k,10},wSUVR_MENG_table(:,2:end));
    end
          
    % Save stats table
    name_excel=fullfile(main_dir,[PETMR_outdata{k,4} '_' tracer '_VOIs_' atlas_mni_name '_' processing_id '.xls']);
    writecell(PETMR_outdata{k,9},name_excel,'Sheet',subjectcode);
    wname_excel=fullfile(main_dir,[PETMR_outdata{k,4} '_' tracer '_VOIs_' atlas_mni_name '_' processing_id '.xls']);
    writecell(PETMR_outdata{k,9},wname_excel,'Sheet',subjectcode);
    
    % Save seg stats table
    name_excel=fullfile(main_dir,['segs_' PETMR_outdata{k,4} '_' tracer '_VOIs_' atlas_mni_name '_' processing_id '.xls']);
    writecell([[cell(1);ROIname],PETMR_outdata{k,10}],name_excel,'Sheet',subjectcode);
    wname_excel=fullfile(main_dir,['segs_' PETMR_outdata{k,4} '_' tracer '_VOIs_' atlas_mni_name '_' processing_id '.xls']);
    writecell([[cell(1);ROIname],PETMR_outdata{k,10}],wname_excel,'Sheet',subjectcode);
    
    if k==4
        % Save colormaps patient space
        Vrefcolormap.fname=fullfile(subjdir_PETcm,['GM_colormap_' atlas_mni_name '_' processing_id '.nii']);
        spm_write_vol(Vrefcolormap,GM_colormap);
        Vrefcolormap.fname=fullfile(subjdir_PETcm,['WM_colormap_' atlas_mni_name '_' processing_id '.nii']);
        spm_write_vol(Vrefcolormap,WM_colormap);
        Vrefcolormap.fname=fullfile(subjdir_PETcm,['CSF_colormap_' atlas_mni_name '_' processing_id '.nii']);
        spm_write_vol(Vrefcolormap,CSF_colormap);
        if WMHC==1
            Vrefcolormap.fname=fullfile(subjdir_PETcm,['WMH_colormap_' atlas_mni_name '_' processing_id '.nii']);
            spm_write_vol(Vrefcolormap,WMH_colormap);
        end
    elseif k==8
        % Save colormaps mni space
        wVrefcolormap.fname=fullfile(subjdir_PETcm,['wGM_colormap_' atlas_mni_name '_' processing_id '.nii']);
        spm_write_vol(wVrefcolormap,wGM_colormap);
        wVrefcolormap.fname=fullfile(subjdir_PETcm,['wWM_colormap_' atlas_mni_name '_' processing_id '.nii']);
        spm_write_vol(wVrefcolormap,wWM_colormap);
        wVrefcolormap.fname=fullfile(subjdir_PETcm,['wCSF_colormap_' atlas_mni_name '_' processing_id '.nii']);
        spm_write_vol(wVrefcolormap,wCSF_colormap);
        if WMHC==1
            wVrefcolormap.fname=fullfile(subjdir_PETcm,['wWMH_colormap_' atlas_mni_name '_' processing_id '.nii']);
            spm_write_vol(wVrefcolormap,wWMH_colormap);
        end
    end
         
end

%%      16/ Save outcome variables and close logfile
%       ------------------------------------

% Define out_name before subjectcode is cleared
out_name=['outcome_variables_' tracer '_' subjectcode '_' processing_id '.mat'];

% Clear variables we don't need to save
cd(subjdir)
clear(list_settings_variables{:})
clear('*img')

% Save workspace variables
save(out_name)

% Save key workspace variables to a csv file
outcome_data = load(out_name);
key_data = array2table(outcome_data.PETMR_outdata(:,[1:8,11:12]),'VariableNames',{'PET','nr_voxels_refVOI','refVOI_value','image','IQR','mean_WM_value','TIV','mean_thickness','bodyweight','injdose'});
writetable(key_data,['key_' out_name(1:end-4) '.csv']);

% Close log file
fprintf(fid,'\n');
fprintf(fid,'Processing finished');
fclose(fid); % close log file

end
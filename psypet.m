function [SUVR_path, SUVR_table_path, SUV_rr_table_path]=psypet(subj, T1, PET, rr, atlas, outfolder)

%% Script written by Thomas Vande Casteele

% Input arguments:
%   1/ subj : the name of the subj (string)
%   2/ T1 : the path to the raw dicom folder, or path to the processed nifti file
%               * if raw dicom folder : we will process the dcm (conversion
%               to nifti, centering, segmentation)
%               * if nifti file : we consider the T1 as already processed
%               by the psypet pipeline. The atlas input argument will be
%               considered as the segmentation and the rr as the path to
%               the reference region in subject space (number or string)
%   3/ PET :  the path to the raw dicom folder, or path to the processed nifti file
%               * if raw dicom folder : we will process the dcm (conversion
%               to nifti, realignement, averaging, SUV calculation)
%               * if nifti file : we consider the PET acquisition as already processed
%               by the psypet pipeline.
%   4/ rr : reference region, can be a number (5), a numeric array ([5,7,8]) or a path to mni (spm) ref region
%   template
%   5/ atlas/segmentation :
%               * if T1 is a nifti file, we consider "atlas" as the result of the segmentation (i.e. aparc+aseg.nii or other)
%                * if T1 is a dcm folder, we consider "atlas" as the one to be obtained: 'Freesurfer', 'Fastsurfer' or 'any CAT12 suppported atlas'
%   (ex: neuromorphometrics)
%   6/ outfolder

%% 0/ Grab script path
script = mfilename('fullpath');
[script_dir,~,~]=fileparts(script);

%% 1/ Process PET if not done yet
if endsWith(PET,'.nii')
    SUV_path=PET;
else
    
    % Make a subfolder
    cd(outfolder)
    mkdir('PETpreproc')
    outfolder_pet=fullfile(outfolder,'PETpreproc');
    
    % Give a name to the processed image
    SUVname=['SUV_' subj];
    
    % convert dcm to nifti and calculate SUV
    [~,~,~,~,~,~,~,filelist_SUV]=LTNP_dcm2SUV(PET,outfolder_pet,SUVname); 
    
    % realign and summate pet frames
    [SUV_path]=LTNP_spm12_realign_and_summate(filelist_SUV,outfolder_pet);  
    
    % Center image
    LTNP_center(SUV_path)
    
end

%% 2/ Process T1 if not done yet
if endsWith(T1,'.nii')
    RR_ready=true;
    T1_path=T1;
    segmentation=atlas;
else
    RR_ready=false;
    
    % Make a subfolder
    cd(outfolder)
    mkdir('T1preproc')
    outfolder_t1=fullfile(outfolder,'T1preproc');
    
    % Give a name to the (to be) processed image
    T1name=['T1_' subj '.nii'];
    
    % dcm2nii, crop, center
    [accT1nii, ~, ~]=LTNP_preproc_T1_spm(T1,outfolder_t1,T1name); 
    
    % delete qform
    LTNP_delete_qform(accT1nii,accT1nii) 
    
    % segment
    if isequal(atlas,'FreeSurfer')
        FreeS=true;
        FastS=false;
        CAT12=false;
        [T1_path,segmentation]=psyfs(accT1nii,outfolder_t1); 
    elseif isequal(atlas,'FastSurfer')
        [T1_path,segmentation]=psyfas(accT1nii,outfolder_t1);
        FreeS=false;
        FastS=true;
        CAT12=false;
    else
        T1_path=accT1nii;
        [segmentation,invdef]=psycat(accT1nii,atlas,outfolder_t1);
        FreeS=false;
        FastS=false;
        CAT12=true;
    end
end

% Grab voxelsize
voxelsize=LTNP_get_voxelsize(T1_path);

%% 3/ Coregister processed PET to processed T1
[rSUV,~]=LTNP_spm12_coregister_reslice(T1_path,SUV_path,outfolder);

%% 4/ Apply RBV PVC on coregistred processed PET
[rbvSUV]=LTNP_PVC_RBV('SUV',script_dir,rSUV,segmentation,outfolder);
%LTNP_PVC_MG()

%% 5/ Make refVOI
if RR_ready
    if ischar(rr)
        refVOI=rr;
    elseif isnumeric(rr)
        refVOI=fullfile(outfolder,'refVOI_mask.nii');
        LTNP_binarize_atlas(segmentation,refVOI,rr,rr);
    end
else
    if FreeS || FastS
        if ischar(rr)
            [~,invdef]=LTNP_calc_deformation_field(T1_path,outfolder);
            [refVOI]=LTNP_spm12_warp_ROI(rr,invdef,outfolder,voxelsize);
        elseif isnumeric(rr)
            refVOI=fullfile(outfolder,'refVOI_mask.nii');
            LTNP_binarize_atlas(segmentation,refVOI,rr,rr);
        end
    elseif CAT12
        if ischar(rr)
            [refVOI]=LTNP_spm12_warp_ROI(rr,invdef,outfolder,voxelsize);
        elseif isnumeric(rr)
            refVOI=fullfile(outfolder,'refVOI_mask.nii');
            LTNP_binarize_atlas(segmentation,refVOI,rr,rr);
        end
    end 
end
[~,refVOI_name,~]=fileparts(refVOI);

%% 6/ Make SUVR image
[~,~,~,SUVR_path]=LTNP_calculate_SUVR(rSUV,refVOI,outfolder); % before rbv
[~,SUVR_name,~]=fileparts(SUVR_path);
[~,~,~,SUVR_rbv_path]=LTNP_calculate_SUVR(rbvSUV,refVOI,outfolder); % after rbv
[~,SUVR_rbv_name,~]=fileparts(SUVR_rbv_path);

%% 7/ Get PET and T1 statistics
[SUVR_table]=LTNP_VOI_stats_v8(SUVR_path,segmentation,'');
[SUV_rr_table]=LTNP_VOI_stats_v8(rSUV,refVOI,'');
[SUVR_rbv_table]=LTNP_VOI_stats_v8(SUVR_rbv_path,segmentation,'');
[SUV_rbv_rr_table]=LTNP_VOI_stats_v8(rbvSUV,refVOI,'');

% Save statistics
SUVR_table_path=fullfile(outfolder,[SUVR_name '.xlsx']);
SUV_rr_table_path=fullfile(outfolder,[refVOI_name '.xlsx']);
SUVR_rbv_table_path=fullfile(outfolder,[SUVR_rbv_name '.xlsx']);
SUV_rbv_rr_table_path=fullfile(outfolder,['PVC_' refVOI_name '.xlsx']);

writecell(SUVR_table,SUVR_table_path,'sheet',subj)
writecell(SUV_rr_table,SUV_rr_table_path,'sheet',subj)
writecell(SUVR_rbv_table,SUVR_rbv_table_path,'sheet',subj)
writecell(SUV_rbv_rr_table,SUV_rbv_rr_table_path,'sheet',subj)

end


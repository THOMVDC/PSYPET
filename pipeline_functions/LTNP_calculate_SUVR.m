function [SUVR, nr_voxels_refVOI, ref_value, SUVR_path]=LTNP_calculate_SUVR(PET,refVOI_thresholded,outfolder,brain_mask_thresholded)

% REFvoi, WMimg, GMimg need to be thresholded (so binary), for example:
%
%   REFvoi_thresholded=refVOI_img > refVOI_thr (0,5)
%   GMimg_thresholded=GMimg > GM_thr (0,3)
%   WMimg_thresholded=sWMimg > WM_thr
%   BRAINmask = (GM_mni + WM_mni + CSF_mni> BRAINmask_thr) (0,95)

% Check input variables
if nargin < 4
    brain_mask_thresholded=1;
end
if isfile(PET)
    [PETimg,Vref]=LCN12_read_image(PET);
else
    PETimg=PET;
end
if isfile(refVOI_thresholded)
    [~, refVOI_name, ~]=fileparts(refVOI_thresholded);
    [refVOI_thresholded]=LCN12_read_image(refVOI_thresholded,Vref);
end

% Apply brain mask to PET data
PETimg=PETimg.*brain_mask_thresholded;

% Create ref_mask
ref_mask = refVOI_thresholded.*brain_mask_thresholded.*(PETimg > 0);

% Get number of voxels and mean intensity in refVOI
nr_voxels_refVOI = sum(ref_mask(:)>0);
ref_value        = nanmean(PETimg(ref_mask>0));

% Calculate SUVR
SUVR             = PETimg./ref_value;

% Save SUVR
if isfile(PET)
    [~, SUV_name, ~]=fileparts(PET);
    SUVR_name=['SUVR_' SUV_name];
    SUVR_path = fullfile(outfolder,[SUVR_name '_' refVOI_name '.nii']); 
    LCN12_write_image(SUVR,SUVR_path,'SUVR',Vref.dt(1),Vref);
else
    SUVR_path='';
end
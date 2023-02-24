function [SUVR, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_pib(PET,refVOI_thresholded,GMimg,GM_thr,brainmask_thresholded)

% REFvoi, WMimg, GMimg need to be thresholded (so binary), for example:
%
%   REFvoi_thresholded=refVOI_img > refVOI_thr (0,5)
%   GMimg_thresholded=GMimg > GM_thr (0,3)
%   WMimg_thresholded=sWMimg > WM_thr
%   BRAINmask = (GM_mni + WM_mni + CSF_mni> BRAINmask_thr) (0,95)

% Apply brain mask to PET data
PET=PET.*brainmask_thresholded;

% Create ref_mask
GMimg_thresholded=GMimg > GM_thr;
ref_mask = refVOI_thresholded.*brainmask_thresholded.*GMimg_thresholded;

% Get number of voxels and mean intensity in refVOI
nr_voxels_refVOI = sum(ref_mask(:)>0);
ref_value        = nanmean(PET(ref_mask>0));

% Calculate SUVR
SUVR             = PET./ref_value;

end

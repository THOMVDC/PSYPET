function [SUVR, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_ucbj(PET,refVOI_thresholded,WMimg,WM_thr,brainmask_thresholded, WMHimg)

% REFvoi, WMimg, GMimg need to be thresholded (so binary), for example:
%
%   REFvoi_thresholded=refVOI_img > refVOI_thr (0,5)
%   GMimg_thresholded=GMimg > GM_thr (0,3)
%   WMimg_thresholded=sWMimg > WM_thr
%   BRAINmask = (GM_mni + WM_mni + CSF_mni> BRAINmask_thr) (0,95)

% Defaults
WMH_thr=0.5;

% Invert WMH image
if isempty(WMHimg)
    iWMHimg=ones(size(WMimg));
else
    WMHimg=WMHimg>WMH_thr;
    iWMHimg=1-WMHimg;
end

% Get the WMH out of the WM
WMimg=WMimg.*iWMHimg;

% Apply brain mask to PET data
PET=PET.*brainmask_thresholded;

% Create ref_mask
sWMimg=zeros(size(WMimg));
spm_smooth(WMimg,sWMimg,[7 7 7]); 
WM_ref_img = sWMimg.*refVOI_thresholded;
WM_ref_img_delta = 0.99*(max(WM_ref_img(:))-min(min(min((WM_ref_img(WM_ref_img(:)>0)))))); % 99% 
WM_ref_mask = WM_ref_img > WM_ref_img_delta;
WMimg_thresholded=WMimg > WM_thr;
ref_mask = (WMimg_thresholded).*WM_ref_mask;     

% Get number of voxels and mean intensity in refVOI
nr_voxels_refVOI = sum(ref_mask(:)>0);
ref_value        = nanmean(PET(ref_mask>0));

% Calculate SUVR
SUVR             = PET./ref_value;

end
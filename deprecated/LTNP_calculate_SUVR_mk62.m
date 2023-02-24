function [SUVR, ref_mask, nr_voxels_refVOI, ref_value]=LTNP_calculate_SUVR_mk62(PET,refVOI_thresholded,GMimg, GM_thr,brainmask_thresholded)

% REFvoi, WMimg, GMimg need to be thresholded (so binary), for example:
%
%   REFvoi_thresholded=refVOI_img > refVOI_thr (0,5)
%   GMimg_thresholded=GMimg > GM_thr (0,3)
%   WMimg_thresholded=sWMimg > WM_thr
%   BRAINmask = (GM_mni + WM_mni + CSF_mni> BRAINmask_thr) (0,95)


% Apply brain mask to PET data
PET=PET.*brainmask_thresholded;

% Smooth
stPET=zeros(size(PET));
spm_smooth(PET,stPET,[4 4 4]);

% Create ref_mask
GMimg_thresholded=GMimg > GM_thr;
ref_mask = brainmask_thresholded.*refVOI_thresholded.*(GMimg_thresholded);

% Get number of voxels and mean intensity in refVOI
nr_voxels_refVOI = sum(ref_mask(:)>0);
ref_value        = nanmean(PET(ref_mask>0));

% Calculate SUVR
SUVR             = PET./ref_value;

end


%BRAINmask_participant=GMimg > GM_threshold;
%sGMimg=zeros(size(GMimg));
%spm_smooth(GMimg,sGMimg,[5 5 5]); 
%GM_ref_img = sGMimg.*refVOIimg;
%GM_ref_img = BRAINmask.*refVOIimg.*(GMimg>GM_threshold);
%GM_ref_img_delta = 0.50*(max(GM_ref_img(:))-min(min(min((GM_ref_img(GM_ref_img(:)>0)))))); % 50% 
%GM_ref_mask = GM_ref_img > GM_ref_img_delta;
%ref_mask = BRAINmask_participant.*GM_ref_mask;
% Smooth MK
%CSF_mask = BRAINmask.*(CSFimg>0.1);
%GM_WM_mask = (GMimg + WMimg) > 0.2;
%BRAINmask_participant=(CSF_mask + GM_WM_mask)>0.9;
%[subjdir_PET,PET_name,PET_ext]=fileparts(SUV);
% stPET_name=['stw' PET_name PET_ext];
% stPET_path=fullfile(subjdir_PET,[stPET_name '.nii']);
% LCN12_write_image(stPET,stPET_path,['masked smoothed' tracer],Vref.dt(1),Vref);
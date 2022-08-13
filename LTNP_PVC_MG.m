function [pvc_MGorig_path,pvc_MGmodif_path,mean_WM_value] = LTNP_PVC_MG(PET_path,GM_path,WM_path,FWHM,out_folder,WM_VOI)

% Background
% ----------
% PVC postprocessing based on the method of Muller-Gartner
%
% input:
%  image = 3D matrix representing a PET/SPECT image.
%  voxelsize = vector of voxelsizes in each direction in mm
%  GM    = 3D matrix representing the corresponding GM fuzzy segmentation
%  WM    = 3D matrix representing the corresponding WM fuzzy segmentation
%  FWHM  = FHWM of the image resolution (assumed isotropic) in mm
%  WM_VOI = 3D matrix of WM VOI to be used. if not specified or empty, we
%           will calculate the WM VOI by taking all WM voxels above a 
%           certain threshold and then erode this VOI. (optional argument) 
% output: 
%  I_PVC = 3D matrix of the partial volume corrected PET/SPECT image using
%          an improved calculation of MG
%  I_PVC_orig = 3D matrix of the partial volume corrected PET/SPECT image 
%               using the original implementation of MG
%  I_WM = estimated intensity in the WM

% Defining parameters
% --------------------
th1      = 0.50; % only voxels with at least this value in the segmentation maps will be analyzed.
th3      = 0.90; % this defines the initial set of voxels in the WM to determine the average true WM value.
min_number_voxels = 20; % the minimum number of voxels to determine the WM activity.
min_th3 = 0.5; % the th3 threshold will not be lowered automatically below this value
kernel  = [3 3 3]; % 5 voxels in each dimension for the erosion step

% Read data
% ---------
% Read PET, get Vref and voxelsize
[PETimg, Vref] =LCN12_read_image(PET_path);
tmp = spm_imatrix(Vref.mat);   % tmp(7:9) are the voxel sizes
voxelsize = abs(tmp(7:9)); % get voxelsize 

% Read GM,WM images with Vref from PET
GMimg=LCN12_read_image(GM_path,Vref);
WMimg=LCN12_read_image(WM_path,Vref);


% Partial Volume Correction
% ---------------------------

% initialisations
meanPET_PVC_orig = zeros(size(PETimg));
meanPET_PVC_modif= zeros(size(PETimg));
cGMimg = zeros(size(PETimg));
cWMimg = zeros(size(PETimg));

% calculate GM mask 
GM_mask  = (GMimg > th1);

% calculate WM mask by excluding WM lesions
if nargin < 6
   WM_VOI = [];
end
if isempty(WM_VOI)
    % calculate 'pure' WM region far enough from GM (for which we determine
    % the WM value estimate). 
    WM_tmp = (WMimg > th3); % initial WM roi to determine average PET value
    WM_VOI = LCN_3Dimage_erode(WM_tmp,kernel); % we erode WM_tmp to avoid to be to close to GM but at the same time it should contain enough voxels. 
    % Check if there are at least min_number_voxels in WM_VOI. If not, use a more liberal threshold and redo.
    nr_voxels = sum(WM_VOI(:));
    while ((nr_voxels < min_number_voxels) && (th3 > min_th3)) % bug fix: && instead of ||
        th3 = 0.9*th3; 
        WM_tmp = (WMimg > th3); % initial WM roi to determine average PET value
        WM_VOI = LCN_3Dimage_erode(WM_tmp,kernel);
        nr_voxels = sum(WM_VOI(:));
    end
end

% determine average WM value in the PET (or SPECT) image
mean_WM_value  = mean(mean(mean(PETimg(WM_VOI))));

% convolution with PSF
smooth_kernel = [FWHM/voxelsize(1) FWHM/voxelsize(2) FWHM/voxelsize(3)];
spm_smooth(1.*GMimg,cGMimg,smooth_kernel); 
spm_smooth(1.*WMimg,cWMimg,smooth_kernel); 

% calculate temporary image (corrects for spill in)
meanPET_tmp = PETimg - cWMimg.*mean_WM_value; 

% original implementation for the calculation of PVC
meanPET_PVC_orig(GM_mask) = meanPET_tmp(GM_mask)./cGMimg(GM_mask);

% improved implementation for the calculation of PVC after deconvolution
sigma = mean(smooth_kernel)./(2.*sqrt(2.*log(2)));
Gauss3Dfilter = fspecial3('gaussian',round(6*sigma),sigma);
J = deconvreg(meanPET_tmp,Gauss3Dfilter); 
meanPET_PVC_modif(GM_mask) = J(GM_mask)./GMimg(GM_mask);
 
% there might be some voxels which actually have negative values in the I_PVC
neg_index = (meanPET_PVC_modif < 0);
meanPET_PVC_modif(neg_index) = 0;

% Save data
[~, PET_name, PET_ext]=fileparts(PET_path);
pvc_MGorig_path = fullfile(out_folder,['PVC_MGorig_' PET_name PET_ext]); 
pvc_MGmodif_path = fullfile(out_folder,['PVC_MGmodif_' PET_name PET_ext]); 
LCN12_write_image(meanPET_PVC_orig,pvc_MGorig_path,'PVC_orig',Vref.dt(1),Vref); 
LCN12_write_image(meanPET_PVC_modif,pvc_MGmodif_path,'PVC_modif',Vref.dt(1),Vref); 

end

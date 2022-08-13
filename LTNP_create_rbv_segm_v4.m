function [output_segm_and_skull_path,MENGSKULL_path,GMWMCSFmask_path,output_segm_path]=LTNP_create_rbv_segm_v4(GMmask_path, WMmask_path, CSFmask_path, atlas_path, output_folder)

% Function written for cat12 segmentations
    % p0 = GM+WM+CSF
    % p1 = GM
    % p2 = WM
    % p3 = CSF
    % p4 = meninges
    % p5 = skull
    % p6 = background 
    
% Load masks and atlas
GM  = LCN12_read_image(GMmask_path);
WM  = LCN12_read_image(WMmask_path);
CSF = LCN12_read_image(CSFmask_path);
[ATLASimg,Vref]   = LCN12_read_image(atlas_path);
ATLASimg = round(ATLASimg);

% Make an initial BRAINmask
BRAINmask=1.*((GM+WM+CSF)>0);

% Smooth brain mask
SBRAINmask=zeros(size(BRAINmask));
spm_smooth(BRAINmask,SBRAINmask,[5 5 5]./[sqrt(8*log(2)) sqrt(8*log(2)) sqrt(8*log(2))],0);

% Make mengeal-skull mask
SBRAINmask=1.*(SBRAINmask>0);
MSmask=SBRAINmask-BRAINmask;

% Create binary masks (logical arrays)
%[GMmask,~,WMmask,~,CSFmask,~,~,~]=LTNP_make_labelimage(GM_path,WM_path,CSF_path);

% Grab voi_values
voi_values = unique(ATLASimg(ATLASimg>0)); % 

% Create ROIs
ATLAS_GM=GM.*ATLASimg;
ATLAS_WM=(voi_values(end)+1)*WM;
ATLAS_CSF=(voi_values(end)+2)*CSF;
ATLAS_MENSKULL=(voi_values(end)+3)*MSmask;

% Include GM rests into ATLAS_CSF
%ATLAS_GM_rest=BRAINmask.*((GM-(ATLASimg>0))>0); % GM not overlapping with an atlas region will be considered as CSF.
%ATLAS_CSF=ATLAS_CSF+((voi_values(end)+2)*(ATLAS_GM_rest));

% Add the new ROIs on the atlas
SEGMENTATION=ATLAS_GM + ATLAS_WM + ATLAS_CSF;
SEGMENTATION_AND_SKULL= ATLAS_GM + ATLAS_WM + ATLAS_CSF + ATLAS_MENSKULL;


% Define paths to save to
[~,atlas_name,~]=fileparts(atlas_path);
output_segm_path=fullfile(output_folder,['rbv_segm_' atlas_name '.nii']);
output_segm_and_skull_path=fullfile(output_folder,['rbv_segm_skull_' atlas_name '.nii']);
MENGSKULL_path=fullfile(output_folder,['mengskull_' atlas_name '.nii']);
GMWMCSFmask_path=fullfile(output_folder,['GMWMCSFmask_' atlas_name '.nii']);

% Save
Vref.fname=output_segm_path;
spm_write_vol(Vref,SEGMENTATION); 
Vref.fname=output_segm_and_skull_path;
spm_write_vol(Vref,SEGMENTATION_AND_SKULL); 
Vref.fname=MENGSKULL_path;
spm_write_vol(Vref,MSmask);
Vref.fname=GMWMCSFmask_path;
spm_write_vol(Vref,BRAINmask); 

end

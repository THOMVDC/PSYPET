function refVOI=make_SO_with_fs(T1_fs, T1, WM, ref_VOI_pt, WML, FLAIR, output_folder)


%% 1/ Coreg FLAIR and WML to T1_fs
if ~isempty(FLAIR)
    [~,rWML]=LTNP_spm12_coregister_reslice(T1_fs,FLAIR,output_folder,WML);
    manWML='';
else
    manWML=WML;
end

%% 2/ Coreg T1,WM, refVOI to T1_fs
other_images={WM;ref_VOI_pt;manWML};
[~,rother_images]=LTNP_spm12_coregister_reslice(T1_fs,T1,other_images);
rWM=rother_images{1};
rref_VOI=rother_images{2};
rmanWML=rother_images{3};
if isempty(FLAIR)
    rWML=rmanWML;
end

%% 3/ Calculate ref_VOI_fs

% Read images
[WMimg,Vref]=LCN12_read_image(rWM);
WMHimg=LCN12_read_image(rWML,Vref);
refVOIimg=LCN12_read_image(rref_VOI);

% Defaults
WMH_thr=0.5;
ref_thr=0.5;

% Threshold ref_VOI
refVOI_thresholded=refVOIimg>ref_thr;

% Invert WMH image
WMHimg=WMHimg>WMH_thr;
iWMHimg=1-WMHimg;

% Get the WMH out of the WM
WMimg=WMimg.*iWMHimg;

% Create ref_mask
sWMimg=zeros(size(WMimg));
spm_smooth(WMimg,sWMimg,[7 7 7]); 
WM_ref_img = sWMimg.*refVOI_thresholded;
WM_ref_img_delta = 0.99*(max(WM_ref_img(:))-min(min(min((WM_ref_img(WM_ref_img(:)>0)))))); % 99% 
WM_ref_mask = WM_ref_img > WM_ref_img_delta;
WMimg_thresholded=WMimg > WM_thr;
ref_mask = (WMimg_thresholded).*WM_ref_mask;     
refVOI=fullfile(output_folder,'ref_name');

% Write image
LCN12_write_image(ref_mask,ref_mask_path,Vref,Vref.dt(1));

end
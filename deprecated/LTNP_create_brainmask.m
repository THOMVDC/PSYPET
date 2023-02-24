function [MASK_path]=LTNP_create_brainmask(GM_path,WM_path,CSF_path,out_folder)

% Load probability maps
[GM,Vref]  = LCN12_read_image(GM_path);
WM         = LCN12_read_image(WM_path);
CSF        = LCN12_read_image(CSF_path);

% Make an initial GMWMCSF mask
MASK=1*((GM+WM+CSF)>0);

% Save mask
[~,name,~]=fileparts(GM_path);
MASK_path=fullfile(out_folder,['GMWMCSF_mask_' name '.nii']);
Vref.fname=MASK_path;
spm_write_vol(Vref,MASK);

end

%BRAINmask_path = (GM_img + WM_img + CSF_img> 0.95) ;
%BRAINmask_mni = (wGM_img + wWM_img + wCSF_img> 0.95) ;

% BET patient (brain extraction tool from FSL) 
% betFILE=fullfile(subjdir_ANAT,'BET',['bet' ref_image_name ref_image_ext]);
% tempfolder=fullfile(subjdir_ANAT,'BET');
% unix(['source ~/.bashrc && cd ' tempfolder ' && bet ' accT1nii ' ' betFILE ' -R -m -n -o -s']);
% mask=fullfile(subjdir_ANAT,'BET',['bet' ref_image_name '_mask' ref_image_ext '.gz']);
% gunzip(mask)
% BRAINmask_path=mask(1:end-3);
% 
% % BET mni
% inFILE=fullfile(subjdir_ANATsegm,'mri',['wm' ref_image_name ref_image_ext]);
% betFILE=fullfile(subjdir_ANAT,'BET',['bet' ref_image_name ref_image_ext]);
% system(['source ~/.bashrc && cd ' tempfolder ' && bet ' inFILE ' ' betFILE ' -R -m -n']);
% mask=fullfile(subjdir_ANAT,'BET',['bet' ref_image_name '_mask' ref_image_ext '.gz']);
% gunzip(mask)
% BRAINmask_mni=mask(1:end-3);
        
% %  HD-BET
%    system(['source ~/.bashrc && cd BET && hd-bet -i ' inFILE ' -o BET -device cpu -mode fast -tta 0']);
%    maskbstFILE=[inFILE '_bet_mask.nii.gz'];
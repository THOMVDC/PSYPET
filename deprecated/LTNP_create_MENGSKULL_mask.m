function [MENGSKULL_path,mask_path]=LTNP_create_MENGSKULL_mask(GM_path,WM_path,CSF_path,out_folder)

% Create brain mask
[mask_path]=LTNP_create_brainmask(GM_path,WM_path,CSF_path,out_folder);

% Smooth
smask=fullfile(out_folder,'stmp_mask.nii');
spm_smooth(mask_path,smask,[5 5 5]./[sqrt(8*log(2)) sqrt(8*log(2)) sqrt(8*log(2))],0);

% Read mask & smask
[MASK,Vref]=LCN12_read_image(mask_path);
SMASK=LCN12_read_image(smask,Vref);

% Get voxels
SMASK=1*(SMASK>0);
MASK=SMASK-MASK;

% Save mask
Vref.fname=fullfile(out_folder,'mengskull.nii');
spm_write_vol(Vref,MASK);

% Return path to mengskull
MENGSKULL_path=Vref.fname;

%% By inflation

% % Initiation
% n=1;
% test=MASK.*BCKG;
% 
% % Dilate mask untill mask (MASK) overlaps with background (BCKG)
% while sum(test(:)) == 0
%     n=n+1;
%     mask_out=fullfile(main_dir,['maskinfl_' num2str(n) '.nii']);
%     unix(['source ~/.bashrc && maskfilter ' mask_in ' dilate ' mask_out]);
%     MASK=LCN12_read_image(mask_out);
%     test=MASK.*BCKG;
%     mask_in=mask_out;
% end
% 
% % Delete all masks except n-1
% delete(mask_out)
% for i=1:(n-2)
%     delete(fullfile(main_dir,['maskinfl_' num2str(i) '.nii']));
% end
% 
% % Read mask n-1
% MENGSKULL=LCN12_read_image(fullfile(main_dir,['maskinfl_' num2str(n-1) '.nii']));
% 
% % Substract
% MASK=1*(GWC>0);
% MENGSKULL=MENGSKULL-MASK;
% 
% % Save
% MENGSKULL_out=fullfile(main_dir,'mengskull_mask.nii');
% LCN12_write_image(MENGSKULL, MENGSKULL_out,'meningeal/skull mask by maskinflation');

end

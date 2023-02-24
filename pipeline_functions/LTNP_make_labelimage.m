function [GMmask_path, WMmask_path, CSFmask_path,BRAINmask_path,p0_path]=LTNP_make_labelimage(GM_path,WM_path,CSF_path,out_folder)
                                                 
% Read segmentations
[GM,Vref] = LCN12_read_image(GM_path);
WM = LCN12_read_image(WM_path,Vref);
CSF = LCN12_read_image(CSF_path,Vref);

% Create binary masks (logical arrays)
%Background=1*((GM+WM+CSF)==0);
%GMmask=(GM>=WM&GM>=CSF)-Background; %GM will be at least 1/3;
%WMmask=(WM>GM&WM>=CSF); %WM will be at least 1/3;
%CSFmask=(CSF>GM&CSF>WM); %CSF will be at least 1/3;

% Create binary masks, Laura's way:
% GMmask=GM>0.3;
% CSF=CSF>0.5;
% WM=WM>0;
% CSFmask=CSF-(CSF.*GM); % Csf that doesn't overlap with GM
% WMmask=WM-GM-CSF;

% Create binary masks, Laura's way:
GMmask=GM>0.3;
CSFmask=CSF>0.5;
WMmask=WM>0;
CSFmask=CSFmask-(CSFmask.*GMmask); % Csf that doesn't overlap with GM
WMmask=WMmask-GMmask-CSFmask;

% Create brainmask
BRAINmask=1*(GMmask+WMmask+CSFmask);

% Create differential segments
G=2.*GMmask;
W=3.*WMmask;
C=1.*CSFmask;

% Sum
p0=G+W+C;

% Grab name
[~,GMname,GMext]=fileparts(GM_path);
[~,WMname,WMext]=fileparts(WM_path);
[~,CSFname,CSFext]=fileparts(CSF_path);
Pname=GMname(GMname==WMname);
Pext=GMext;

% Set path
GMmask_path=fullfile(out_folder,['mask_' GMname GMext]);
WMmask_path=fullfile(out_folder,['mask_' WMname WMext]);
CSFmask_path=fullfile(out_folder,['mask_' CSFname CSFext]);
BRAINmask_path=fullfile(out_folder,['GMWMCSF_mask_' Pname Pext]);
p0_path=fullfile(out_folder,['label_' Pname Pext]);

% Save
LCN12_write_image(GMmask,GMmask_path,'GMmask',Vref.dt(1),Vref);
LCN12_write_image(WMmask,WMmask_path,'WMmask',Vref.dt(1),Vref);
LCN12_write_image(CSFmask,CSFmask_path,'CSFmask',Vref.dt(1),Vref);
LCN12_write_image(BRAINmask,BRAINmask_path,'BRAINmask',Vref.dt(1),Vref);
LCN12_write_image(p0,p0_path,'labelimage',Vref.dt(1),Vref);

end
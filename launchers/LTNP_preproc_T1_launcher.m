
%% Background

% Launcher example for LTNP_preproc_T1_spm
%
% Author: 
%       Thomas Vande Casteele, KU Leuven

%% Settings

% Define subject, input directory, output directory and the output prefix
subj='B001';
inputdir = '/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/UCBJ_4s/';
outputdir = '/Users/tvdcas1/Desktop/';
outputprefix = 'T1';


%% Launch for single subject

% Assemble dicom directory
dcmdir=fullfile(inputdir,subj);

% Assemble output name
outputname=fullfile(outputdir,[outputprefix '_' subj '.nii']);

% Launch
[~,~,~]=LTNP_preproc_T1_spm(dcmdir,outputdir,outputname); 


%% Launch for multiple subjects

% List subjects
l=dir(inputdir);

for s=1:length(subjects)
    
    % Grab subject
    subj=l(s).name;
    
    % Assemble dicom directory
    dcmdir=fullfile(inputdir,subj);
    
    % Assemble outputname
    outputname=fullfile(outputdir,[outputname '_' subj '.nii']);
    
    % Launch
    [~,~,~]=LTNP_preproc_T1_spm(dcmdir,outputdir,outputname); 
end



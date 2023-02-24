
%% Background

% Launcher example for LTNP_dcm2SUV 
%
% Author: 
%       Thomas Vande Casteele, KU Leuven

%% Settings

% Define subject, input directory, output directory and the output prefix
subj='B001';
inputdir = '/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET/UCBJ_4s/';
outputdir = '/Users/tvdcas1/Desktop/';
outputprefix = 'SUV';


%% Launch for single subject

% Define dicom directory
dcmdir=fullfile(inputdir,subj);

% Assemble dicom directory
outname=fullfile(outputdir,[outputprefix '_' subj]);

% Launch
[~,~]=LTNP_dcm2SUV(dcmdir,outputdir,outname);


%% Launch for multiple subjects

% List subjects
l=dir(inputdir);

for s=1:length(l)
    
    % Grab subject
    subj=l(s).name;
    
    % Assemble dicom directory
    dcmdir=fullfile(inputdir,subj);
    
    % Assemble output name
    outname=fullfile(outputdir,[outputname '_' subj]);
    
    % Launch
    [~,~]=LTNP_dcm2SUV(dcmdir,outputdir,outname);
end


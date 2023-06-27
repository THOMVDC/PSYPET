%% Background

% Launcher example for LTNP_preproc_T1_spm
%
% Author: 
%       Thomas Vande Casteele, KU Leuven


%% Settings

% Define subject, input directory and the output path for excel
subj='B002';
inputdir='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/MK62_4s';
output_excel_path = '/Volumes/LaCie/Thomas/Projects/MK_paper/clinical data/SUV_info_update.xlsx';

%% Launch for single subject

% Assemble dicom directory
dcmdir=fullfile(inputdir,subj);

% Define table 
T = cell(1,9);
T{1,1} = 'subject';
T{2,1} = subj;
T{1,1+1}= 'ptname';
T{1,1+2}= 'ptid';
T{1,1+3} =  'bodyWeight';
T{1,1+4} =  'injdosis';
T{1,1+5} =  'acqtime';
T{1,1+6} =  'injtime';
T{1,1+7} =  'halftime';
T{1,1+8} =  'rescale_factor';

% Launch
[T{2,2},T{2,3},T{2,4},T{2,5},T{2,6},T{2,7},T{2,8},T{2,9}]=LTNP_SUV_info_retriever(dcmdir);

% Write subject table to excel
writecell(T,output_excel_path);


%% Launch for multiple subjects

% List all subjects
l = dir(inputdir);
nr_subj=length(l);

% Define table 
T = cell(nr_subj,9);
T{1,1} = 'subject';
T{1,1+1}= 'ptname';
T{1,1+2}= 'ptid';
T{1,1+3} =  'bodyWeight';
T{1,1+4} =  'injdosis';
T{1,1+5} =  'acqtime';
T{1,1+6} =  'injtime';
T{1,1+7} =  'halftime';
T{1,1+8} =  'rescale_factor';

for i = 1:nr_subj
    
    % Grab subject
    subj = l(i).name; 
    T{1+i,1} = subj;
    
    % Assemble dicom directory
    dcmdir=fullfile(l(i).folder,subj);
    
    % Launch
    [T{1+i,2},T{1+i,3},T{1+i,4},T{1+i,5},T{1+i,6},T{1+i,7},T{1+i,8},T{1+i,9}]=LTNP_SUV_info_retriever(dcmdir);
    
end

% Write the population table to excel
writecell(T,output_excel_path);

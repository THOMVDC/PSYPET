    
% Define in- and output folders
T1dir='/Volumes/LaCie/Thomas/Projects/RETINAL_IMAGING/DATA/RAW/T1';
outfolder='/Volumes/LaCie/Thomas/Projects/RETINAL_IMAGING/DATA/FREESURFER7/preproc';
fsfolder='/Volumes/LaCie/Thomas/Projects/RETINAL_IMAGING/DATA/FREESURFER7/input';

% List subjects
subjects=dir(fullfile(T1dir, '*'));

%% Preproc T1
for s=1:length(subjects)
    
    % Grab subject code
    subj=subjects(s).name;
    
    % Define and make outfolder subject
    outfolder_subj=fullfile(outfolder,subj);
    cd(outfolder)
    mkdir(subj)
    
    % Grab T1
    T1=fullfile(T1dir,subj);
    
    % Give a name to the (to be) processed image
    T1name=['T1_' subj '.nii'];
    
    % dcm2nii, crop, center
    [accT1nii, ~, ~]=LTNP_preproc_T1_spm(T1,outfolder_subj,T1name); 
    
    % delete qform
    %LTNP_delete_qform(accT1nii,accT1nii) 
    
    copyfile(accT1nii,fsfolder)
    
end
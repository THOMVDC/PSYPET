% note à moi meme est-ce qu'en utilisant dcm converter de spm, j'ai pas une
% autre spatialisation que dcm2nii ? c'est un des truc que j'ai changé dans
% LTNP_preproc_T1 qui peut expliquer le résultat différent dans le
% outfolder ci bas

subjects={'B069','B070','B071','B072', 'B073', 'B074', 'B075', 'B076', 'B077'};

for s=1:length(subjects)

    infolder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/T1/DICOM';
    outfolder='/Volumes/LaCie/L3D/L3D_shared/processed/T1_DATA/CAT12';
    subjectcode=subjects{s};
    dcmfolder=fullfile(infolder,subjectcode);
    niifolder=fullfile(outfolder,subjectcode);
    outname=['T1_' subjectcode '.nii'];
    mni_origin='/Volumes/LaCie/Thomas/Projects/SCRIPTS/PSYPET/templates/ACsphere_bin.nii.gz';
    mni='/Volumes/LaCie/Thomas/Projects/SCRIPTS/PSYPET/templates/Template_T1_IXI555_MNI152_GS.nii';
    cd(outfolder)
    mkdir(subjectcode)
    LTNP_preproc_T1(dcmfolder,niifolder,outname,mni,mni_origin)  
end

 

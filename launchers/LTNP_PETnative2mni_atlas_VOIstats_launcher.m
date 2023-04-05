
subjects={'B003','B004','B005','B006','B008','B009','B010','B011','B013','B014','B015','B016','B017','B019','B020','B021','B023','B024','B026','B027','B028','B030','B031','B033','B034','B035','B037','B038','B040','B042','B043','B045','B046','B047','B048','B049','B050','B051','B052','B053','B055','B056','B057','B059','B061','B063','B065','B072'};
SUVR_folder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/UCBJ/output/psypet_cat12_output_SOnoWML';
ATLASmni='/Users/mlaroy0/spm12/toolbox/cat12/templates_1.50mm/neuromorphometrics.nii';
output_folder='/Users/tvdcas1/Desktop/test';
T1_folder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output';

for s=1:length(subjects)
    subj=subjects{s};
    SUVR=fullfile(SUVR_folder,subj,'SUVR_SUV_PET_PVC_RBV_65mm_in_seg_refVOI_without_WML.nii');
    deformation_field_or_T1=fullfile(T1_folder,'cat12_output','mri',['y_accT1_' subj '.nii']);
    output_folder_subj=fullfile(output_folder,subj);
    cd(output_folder)
    mkdir(subj)
    % Make it run !
    [SUVR_mni_table,~]=LTNP_PETnative2mni_atlas_VOIstats(SUVR, ATLASmni, output_folder_subj, deformation_field_or_T1);
end
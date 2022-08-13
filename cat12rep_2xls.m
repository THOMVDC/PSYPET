
maindir_T1='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/cat12_output/report';
out_excel='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/PSYPET4/T1/output/cat12_stats/cat12_L3D_MCI3.xlsx';

% CAT12.7
T_subjects=dir([maindir_T1 '/B0*']);
for s=1:length(T_subjects)
    subj=T_subjects(s).name;
    infile=dir(fullfile(maindir_T1,subj,'psypet*2022','ANAT','CAT12','report',['cat_accT1_' subj '.mat']));
    infile=fullfile(infile.folder,infile.name);
    %infile=fullfile(maindir_T1,subj,'report',['cat_bcaccT1_' subj '.mat']);
    CAT12_vol_thick=load(infile);
    str=CAT12_vol_thick.S.catlog{end-4, 1};
    T_subjects(s).TIV=CAT12_vol_thick.S.subjectmeasures.vol_TIV;
    T_subjects(s).thickness=cell2mat(CAT12_vol_thick.S.subjectmeasures.dist_thickness);
    T_subjects(s).thickness=T_subjects(s).thickness(1); % grab only mean, not confidence interval
    T_subjects(s).IQR=extractBetween(str,'(IQR):  ','%');
    T_subjects(s).CSF=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(1);
    T_subjects(s).GM=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(2);
    T_subjects(s).WM=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(3);
    T_subjects(s).WML=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(4);
end
   
writetable(struct2table(T_subjects), out_excel)

% CAT12.6
T_subjects=dir(fullfile(maindir_T1,'cat_accT1_*.mat'));
for s=1:length(T_subjects)
    %subj=T_subjects(s).name;
    %infile=dir(fullfile(maindir_T1,subj,'psypet*2022','ANAT','CAT12','report',['cat_accT1_' subj '.mat']));
    %infile=fullfile(infile.folder,infile.name);
    infile=fullfile(T_subjects(s).folder,T_subjects(s).name);
    CAT12_vol_thick=load(infile);
    %str=CAT12_vol_thick.S.catlog{end-4, 1};
    %T_subjects(s).name=T_subjects(s).name(11:14);
    T_subjects(s).name=T_subjects(s).name(11:16);
    T_subjects(s).TIV=CAT12_vol_thick.S.subjectmeasures.vol_TIV;
    %T_subjects(s).thickness=cell2mat(CAT12_vol_thick.S.subjectmeasures.dist_thickness);
    %T_subjects(s).thickness=T_subjects(s).thickness(1); % grab only mean, not confidence interval
    %T_subjects(s).IQR=extractBetween(str,'(IQR):  ','%');
    T_subjects(s).CSF=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(1);
    T_subjects(s).GM=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(2);
    T_subjects(s).WM=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(3);
    T_subjects(s).WML=CAT12_vol_thick.S.subjectmeasures.vol_abs_CGW(4);
end

writetable(struct2table(T_subjects), out_excel)

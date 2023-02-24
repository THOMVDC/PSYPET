function [acIMG, D]=LTNP_set_origin_AC(IMG,mni,mni_origin)

% Grab input folder, file, and extension
[infolder, infile, inext]=fileparts(IMG);

% Define temporary variables used by ANTS
output_s='mnitemplateinpatientspace';
invmat_s=[output_s '0GenericAffine.mat'];
mniACspherePat='mniACspherePat.nii.gz';

% Set name and path of your outputimage
acIMG=fullfile(infolder,['ac' infile inext]);

% Make temporary folder in the input folder
cd(infolder);
mkdir('tmp');
tmpfolder=fullfile(infolder,'tmp');

% Bring mni to T1, use the same matrix to bring the ACsphere in patient space
system(['source ~/.bashrc && cd ' tmpfolder ' && antsRegistrationSyN.sh -f ' IMG ' -m ' mni ' -t a -x NULL -o ' output_s]);
system(['source ~/.bashrc && cd ' tmpfolder ' && WarpImageMultiTransform 3 ' mni_origin ' ' mniACspherePat ' -R ' IMG ' ' invmat_s]);

% Get the origin voxels of the ACsphere in patient space
command_3=['source ~/.bashrc && cd ' tmpfolder ' && fslstats ' mniACspherePat ' -C'];
[~,mrvoxels]=system(command_3);
A=regexp(mrvoxels, '\r?\n', 'split');
B=cell2mat(A(1));
C=split(B);
D=cell2mat([C(1) ' ' C(2) ' -' C(3)]);

% Set origins of image to the one of the ACsphere in patientspace
system(['source ~/.bashrc && cd ' tmpfolder ' && SetOrigin 3 ' IMG ' ' acIMG ' ' D]);

% Remove temporary folder
cd(infolder)
rmdir('tmp','s')

end


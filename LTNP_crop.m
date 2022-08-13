function [cIMG]=LTNP_crop(IMG,outfolder)

% Input should be nifti, output will be nifti

% Define in- and output files
[~, infile, inext]=fileparts(IMG);
aIMG=fullfile(outfolder,['a' infile inext]);
bIMG=fullfile(outfolder,['b' infile inext]);
tempIMG=fullfile(outfolder,['temp' infile inext]);
cIMG=fullfile(outfolder,['c' infile inext]);

% Crop IMG with 3dAutobox (AFNI) to aIMG
system(['source ~/.bashrc && 3dAutobox -input ' IMG ' -prefix ' aIMG ' -npad 10']);

% Set voxels of cropped image to 1, save it to bIMG
[a,Vref]=LCN12_read_image(aIMG);
a=ones(size(a));
LCN12_write_image(a,bIMG,'image of ones',Vref.dt(1),Vref)

% Crop bIMG with robustfov (fsl) to cIMG
system(['source ~/.bashrc && robustfov -i ' bIMG ' -r ' tempIMG ]);

% Read your original image in the cropped volume
gunzip([tempIMG '.gz'])
[~,U]=LCN12_read_image(tempIMG);
P=LCN12_read_image(aIMG,U);
LCN12_write_image(P,cIMG,'the end',U.dt(1),U);

% Delete aIMG
delete(aIMG)
delete(bIMG)
delete(tempIMG)

end
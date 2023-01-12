function [voxelsize,dimension]=LTNP_get_voxelsize_and_dimension(image)

[~,~,ext]=fileparts(image);

if isequal(ext,'.nii')
    info = niftiinfo(image);
    voxelsize=info.PixelDimensions;
    dimension=info.ImageSize;
elseif isequal(ext,'.dcm')
    info = dicominfo(image,'UseDictionaryVR',true);
    voxelsize=[info.PixelSpacing(1) info.PixelSpacing(2) info.SliceThickness];
    dimension=[info.Rows info.Columns info.NumberOfSlices];
elseif isfolder(image)
    list=dir(image);
    list=list(~startsWith({list.name}, '.')); % removing hidden files from the list
    dicomimg=fullfile(image,list(1).name); % grabs first image (pseudorandom)
    info = dicominfo(dicomimg,'UseDictionaryVR',true);
    voxelsize=[info.PixelSpacing(1) info.PixelSpacing(2) info.SliceThickness]; 
    dimension=[info.Rows info.Columns info.NumberOfSlices];
else
    error('Input file specified is not a valid .dcm or a .nii file')
end
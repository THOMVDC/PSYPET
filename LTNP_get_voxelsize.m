function voxelsize=LTNP_get_voxelsize(image)

[~,~,ext]=fileparts(image);

if isequal(ext,'.nii')
    info = niftiinfo(image);
    voxelsize=info.PixelDimensions;
elseif isequal(ext,'.dcm')
    info = dicominfo(image,'UseDictionaryVR',true);
    voxelsize=[info.PixelSpacing(1) info.PixelSpacing(2) info.SliceThickness];
elseif isfolder(image)
    list=dir(image);
    list=list(~startsWith({list.name}, '.')); % removing hidden files from the list
    dicomimg=fullfile(image,list(1).name); % grabs first image (pseudorandom)
    info = dicominfo(dicomimg,'UseDictionaryVR',true);
    voxelsize=[info.PixelSpacing(1) info.PixelSpacing(2) info.SliceThickness]; 
else
    error('Input file specified is not a valid .dcm or a .nii file')
end
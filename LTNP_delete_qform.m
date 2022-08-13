function LTNP_delete_qform(img_path,out_path)

if nargin<2
    out_path=img_path;
end

% Read image from nifti file
img=niftiread(img_path);

% Read header from nifti file
info=niftiinfo(img_path);

% Set sform to 0
info.raw.sform_code=0;

% Write image with new header to output nifti file
niftiwrite(img,out_path,info);

end
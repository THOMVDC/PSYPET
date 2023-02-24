function [bin_img,bin_img_path]=LTNP_binarize_atlas(atlas_path,bin_img_path,downthreshold,upthreshold)
% downthreshold=0.3; upthreshold=Inf;
% if down- = upthreshold; then we match the image to the threshold. 
% You can use an array of numbers to match here. Example : [2,17,9]

% Read image
[ATLASimg,Vref] = LCN12_read_image(atlas_path);
IMG = round(ATLASimg);

% Threshold
if isequal(downthreshold,upthreshold)
    bin_img=1.*(ismember(IMG,downthreshold));
else
    bin_img=1.*(downthreshold<IMG).*(IMG<upthreshold);
end

% Save
Vref.fname=bin_img_path;
spm_write_vol(Vref,bin_img); 

end
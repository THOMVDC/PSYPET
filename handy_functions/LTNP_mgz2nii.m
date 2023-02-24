function [nii]=LTNP_mgz2nii(mgz,out_folder)

% this function requires freesurfer installed correctly and available
% through the PATH

[mgz_folder, mgz_file, ~]=fileparts(mgz);
if nargin<2
    nii=fullfile(mgz_folder,[mgz_file '.nii']);
else
    nii=fullfile(out_folder,[mgz_file '.nii']);
end
    
% Execute mri_convert from freesurfer
system(['source ~/.bashrc && mri_convert -i ' mgz ' -o ' nii]);

end

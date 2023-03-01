function [rbv_path]=LTNP_PVC_RB(subject,rbv_script_dir,pet_path,seg_path,output_folder,fwhm)

% Pet_path and seg_path should be coregistred 

% Relies on the script of Nathalie Mertens
% Dependencies:
% To check which python installation matlab calls:
%   pyenv
% To change the python installation matlab calls (example):
%   pyenv('Version','/Library/Frameworks/Python.framework/Versions/3.6/bin/python3.6')

% On iPsychiater, the directory where python rbv script is located is 
% script_dir='/Volumes/LaCie/Thomas/Projects/L3D/SCRIPTS/NATHALIE_MERTENS/'
% script_name='PVC_RBV_calculation_Nathalie.py'

% Add python script_dir to python path
if count(py.sys.path,rbv_script_dir) == 0
    insert(py.sys.path,int32(0),rbv_script_dir)
end

% Add script_dir to matlab path
addpath(rbv_script_dir)

% Make temporary outputfolder
%output_folder_tmp=fullfile(output_folder,'tmp/');
%cd(output_folder);
%mkdir('tmp');
% if ~endsWith(output_folder,'/')
%     output_folder=[output_folder '/'];
% end

% Run script
rbv_path=py.PVC_RB_calculation_Nathalie.rb_pvc(subject,pet_path,seg_path,output_folder,fwhm);

% Convert python string to matlab string
rbv_path=char(rbv_path); 

% remove temporary folder
%rmdir('tmp')

end
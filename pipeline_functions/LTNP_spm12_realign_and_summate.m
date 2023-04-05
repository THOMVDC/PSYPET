function [out]=LTNP_spm12_realign_and_summate(filelist, output_folder)

% Input:
%       Absolute path to the input_folder, containing the PET frames in nifti format to
%       realign and summate. No other files are contained in the
%       input_folder
%       Absolute path (containing the chosen name to the output_image you want the result to be in
%       
% Output:
%       One image resulting from the realignment and summation of the
%       images found in the input folder
%
% Author: 
%       Thomas Vande Casteele

% Make sure the number of frames is more than one
filelist = filelist(~startsWith({filelist.name}, '.'));
nr_frames = length(filelist);
if nr_frames < 2
    out=filelist;
else

    % Put frames into one variable and copy the original files to the ouput directory 
    Pdy = cell(nr_frames,1);
    for frame = 1:nr_frames
        clear filename src dst
        filename = filelist(frame).name;
        infolder=filelist(frame).folder;
        if filename(1)=='.'
            continue
        end
        src=fullfile(infolder,filename);
        dst=fullfile(output_folder,filename);
        if strcmp(src,dst)==0
            copyfile(src,dst);
        end
        Pdy{frame} = dst;       
    end

    % define filelist of frames to realign
    framelist = cell(nr_frames,1);
    for i = 1:nr_frames
       framelist{i,1} = [char(Pdy{i}) ',1'];
    end        

    % Initialise spm_jobman
    spm_jobman('initcfg')

    % Write batchfile for realignment of frames in filelist, realignment is
    % done to the first image
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.realign.estimate.data = {framelist}; %'
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
    matlabbatch{2}.spm.spatial.realign.write.data = framelist;
    matlabbatch{2}.spm.spatial.realign.write.roptions.which = [0 1]; % % write mean image only, not the separate frames
    matlabbatch{2}.spm.spatial.realign.write.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{2}.spm.spatial.realign.write.roptions.prefix = 'r';

    % Save batchfile
    cd(output_folder)
    save batch_realign matlabbatch

    % Run batchfile
    spm_jobman('run',matlabbatch)

    % Display realignment parameters
    figure(1)
    cd(output_folder)
    filelist = dir('rp_*.txt');   % rp_* is created by SPM after realigning
    if size(filelist,1) ~= 1
       fprintf('ERROR no or more than one realignment file (rp_*.txt) found \n');
       return;
    end
    rp_file = filelist(1).name;
    threshold_totmove=1;
    threshold_translation=3;
    threshold_rotation=3;
    LCN12_analyze_headmovement(rp_file,threshold_totmove,threshold_translation,threshold_rotation); % Calculates and displays momentary displacements as function of time.

    % Rename the mean image (Realign job in SPM does produce mean images (realignment is in two steps: all images are realigned to the first and a mean image is computed to which all images are realigned in a second pass). cfr. line 260: matlabbatch{2}.spm.spatial.realign.write.data = filelist;)
    in = dir('mean*.nii');
    if size(filelist,1) ~= 1
        fprintf('ERROR subject: no or more than one mean image found\n');
        return;
    end
    out=fullfile(output_folder,['r' filename(1:end-10) '.nii']);
    movefile(in(1).name,out);

end
    
end

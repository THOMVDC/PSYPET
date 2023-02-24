function [rescale_factor, filelist_SUV]=LTNP_dcm2SUV(dicomdir,outputdir,outname)

    
    %% Get rescale factor
    [~,~,~,~,~,~,~,rescale_factor]=LTNP_SUV_info_retriever(dicomdir);

   
    %% Run dicom converter from spm
    % Make temporary outputdir
    cd(outputdir);
    mkdir('tmp')
    outputdir_TMP=fullfile(outputdir,'tmp');
    
    % Grab dcm frames
    nr_dcmframes=length(list);
    dcmframes = cell(nr_dcmframes,1);
    for i = 1:nr_dcmframes
       dcmframes{i,1} = fullfile(list(i).folder,list(i).name);
    end        

    % Initialise spm_jobman
    spm_jobman('initcfg')

    % Write job
    matlabbatch{1}.spm.util.import.dicom.data = dcmframes;
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = {outputdir_TMP};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.meta = 1;
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

    % Save batch
    cd(outputdir)
    save batch_dcm2nii matlabbatch
    
    % Run batchfile
    spm_jobman('run',matlabbatch)
       
    %% calculate SUV image
    cd(outputdir_TMP)
    niiframes=dir(fullfile(outputdir_TMP,'*.nii'));
    nr_niiframes=length(niiframes);
    [~,Vref]=LCN12_read_image(niiframes(1).name);
    SUV = zeros(Vref.dim(1),Vref.dim(2),Vref.dim(3),nr_niiframes);
    for i = 1:nr_niiframes
        
        % Grab frame
        tmp = LCN12_read_image(niiframes(i).name,Vref);
        
        % Apply rescale_factor
        tmp = tmp * rescale_factor;
        
        % Build SUV image
        SUV(:,:,:,i) = tmp;
        
        % Save Image
        LCN12_write_image(tmp,fullfile(outputdir,[outname '_0000' num2str(i) '.nii']),'SUV',Vref.dt(1),Vref);
        delete(niiframes(i).name)
        
    end
    delete('*.json')
    
    % Bring the filelist back
    filelist_SUV= dir(fullfile(outputdir,[outname '_0000*.nii']));
    
    % Delete temporary folder 
    cd(outputdir)
    rmdir('tmp','s')
    
end
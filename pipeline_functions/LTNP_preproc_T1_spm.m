function [accT1nii, ptname, ptID]=LTNP_preproc_T1_spm(infolder,outfolder,outname)    % !!! Outname should end with '.nii' as extension !!!
    
    % Convert T1 to nii
    accT1nii=fullfile(outfolder,['acc' outname]); % path to the output nifti file
    
    if endsWith(infolder,'.nii') % assumed already transformed to nifti, cropped and set to ACPC
        
        ptname='unknown';
        ptID='unknown';
        
        % Center T1
        LTNP_center(infolder,accT1nii)
        
        % No crop
        
    else
        [ptname, ptID]=LTNP_dcm2nii(infolder,outfolder,outname(1:end-4)); % transform to nifti
        
        % Center T1
        T1nii=fullfile(outfolder,outname);
        LTNP_center(T1nii,accT1nii)
        
        % No crop
        
    end    
end
function [accT1nii, ptname, ptID]=LTNP_preproc_T1_spm(dcmfolder,outfolder,outname)    % !!! Outname should end with '.nii' as extension !!!
    
    % Convert T1 to nii
    accT1nii=fullfile(outfolder,['acc' outname]); % path to the eventual output nifti file
    if endsWith(dcmfolder,'.nii') % assumed already transformed to nifti, cropped and set to ACPC
        copyfile(dcmfolder,accT1nii);
        ptname='unknown';
        ptID='unknown';
        D='unknown';
    else
        [ptname, ptID]=LTNP_dcm2nii(dcmfolder,outfolder,outname(1:end-4)); % transform to nifti
        % Center T1
        T1nii=fullfile(outfolder,outname);
        accT1nii=fullfile(outfolder,['acc' outname]);
        LTNP_center(T1nii,accT1nii)
    end    
end
function [accT1nii, ptname, ptID, D]=LTNP_preproc_T1(dcmfolder,outfolder,outname,mni,mni_origin)    % Outname should end with '.nii' as extension
    
    % Convert T1 to nii
    accT1nii=fullfile(outfolder,['acc' outname]); % path to the eventual output nifti file
    if endsWith(dcmfolder,'.nii') % assumed already transformed to nifti, cropped and set to ACPC
        copyfile(dcmfolder,accT1nii);
        ptname='unknown';
        ptID='unknown';
        D='unknown';
    else
        % [ptname, ptID]=LTNP_dcm2nii(dcmfolder,outfolder,outname(1:end-4)); % transform to nifti with spm
        [ptname, ptID]=LTNP_dcm2niix(dcmfolder,outfolder,outname(1:end-4)); % transform to nifti with dcm2niix
        

        % Crop T1
        T1nii=fullfile(outfolder,outname);
        [cT1nii]=LTNP_crop(T1nii,outfolder);

        % Set T1 origin to the ACPC
        [accT1nii, D]=LTNP_set_origin_AC(cT1nii,mni,mni_origin);

    end
        
end


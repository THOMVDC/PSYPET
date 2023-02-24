function refVOI=make_SO_with_fs_and_template_laura(ref_VOI_patient, WML, output_folder)

% ref_VOI_pt and WML should be in same space
%% 2/ Calculate ref_VOI_fs

% Read images
[refVOIimg,Vref1]=LCN12_read_image(ref_VOI_patient);
[WMHimg,Vref2]=LCN12_read_image(WML);

if isequal(Vref1.dim,Vref2.dim)
    
    % Defaults
    WMH_thr=0.5;
    ref_thr=0.5;

    % Threshold ref_VOI
    refVOI_thresholded=refVOIimg>ref_thr;

    % Invert WMH image
    WMHimg_thresholded=WMHimg>WMH_thr;
    iWMHimg=1-WMHimg_thresholded;

    % Get the WMH out of the WM
    refVOI_without_WMH=refVOI_thresholded.*iWMHimg;

    % Write image
    refVOI=fullfile(output_folder,'refVOI_without_WML.nii');
    LCN12_write_image(refVOI_without_WMH,refVOI,'refVOI_without_WML',Vref1.dt(1),Vref1);
    %LCN12_write_image(iWMHimg,refVOI,'refVOI_without_WML',Vref1.dt(1),Vref1);
end

end
function [SUVRmni_table,SUVRmni_table_path]=LTNP_PETnative2mni_atlas_VOIstats(SUVR, ATLASmni, output_folder, deformation_field_or_T1)

%% Background
% Assumptions
% ATLAS is in MNI152 space as the template
%
%
% Author:
%       Thomas Vande Casteele, KU Leuven
% 
% Dependency:
%       CAT12, Friedrich  Schiller University  Jena,  Jena,  Germany

%% Processing
% Get voxelsize and dimension of ATLAS
[ATLASmni_vs,ATLASmni_dim]=LTNP_get_voxelsize_and_dimension(ATLASmni);

% Get header info of deformation_field_or_T1
header=niftiinfo(deformation_field_or_T1);

% Check if deformation_field_or_T1 is a deformation_field or a T1
if isequal(header.Description,'Deformation') %deformation_field_or_T1 is a deformation field
    deformation_field=deformation_field_or_T1;
else %deformation_field_or_T1 is a T1
    [deformation_field,~]=LTNP_cat12_calc_deformation_field(deformation_field_or_T1,output_folder,ATLASmni_vs);
end

% Warp SUVR to mni
SUVRmni=LTNP_cat12_warp(SUVR,deformation_field,output_folder);

% Get voxelsize and dimension of SUVRmni
[SUVRmni_vs,SUVRmni_dim]=LTNP_get_voxelsize_and_dimension(SUVRmni);

% Calculate ATLASmni stats of SUVRmni if dimensions and voxelsizes of SUVRmni and ATLAS do correspond
if isequal(SUVRmni_vs,ATLASmni_vs) && isequal (SUVRmni_dim,ATLASmni_dim)
    [SUVRmni_table,~,~]=LTNP_VOI_stats(SUVRmni,ATLASmni,'');
else
    error('SUVR warped to mni has not same voxelsize or dimension as ATLASmni, please check voxelsize and dimension of deformation field')
end

% Save statistics
[~,SUVRmni_name,~]=fileparts(SUVRmni);
[~,ATLASmni_name,~]=fileparts(ATLASmni);
SUVRmni_table_path=fullfile(output_folder,[SUVRmni_name '_' ATLASmni_name '.xlsx']);
writecell(SUVRmni_table,SUVRmni_table_path)

end
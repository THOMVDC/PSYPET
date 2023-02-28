function [table_out,colormap,Vref]=LTNP_VOI_stats(img_path,atlas_path,VOIdetails_path,mask_path)
%% Background %%
%%%%%%%%%%%%%%%%

% To launch this function, I refer to to the LTNP_VOI_stats_v8_launcher.m
%
% Input
%   atlas_path can be a segmentation, parcellation or VOI image
%   img_path and atlas_path should be registred in the same space
%   VOIdetails should be either 
%       empty string (''), in which case the atlas values are used as ROInames 
%       or a csv/excel file, with 3 columns:
%           column 1 = ROIid (can be more than one, ex: [3 4 7 8] for a composite VOI)
%           column 2 = any kind
%           column 3 = ROIname
%   mask_path is an optional argument, he should be registred in the same space as
%       atlas_path and img_path
% 
% Output:
%       table_out=table of statistics (mean, median, std, min, max, vol)
%       within each ROI from the input image
%       colormap= map of the regions considered in VOIdetails 
%                 will be equal to the specified atlas if VOIdetails_path is empty or strictly
%                 equivalent to all the atlas VOIs and if no mask is used
%
% Author: 
%       Thomas Vande Casteele, KU Leuven
%
% Dependency:
%       SPM12, Wellcome Trust Centre for Neuroimaging, University College London


%% Processing %%
%%%%%%%%%%%%%%%%

% Read atlas and img
[ATLASimg,Vref] = LCN12_read_image(atlas_path);
ATLASimg = round(ATLASimg);
IMG= LCN12_read_image(img_path,Vref); 
[~,img_name,~]=fileparts(img_path);

% Read mask
if nargin<4
    MASKimg=1;
    msk_name='';
else
    MASKimg=LCN12_read_image(mask_path,Vref); 
    [~,msk_name,~]=fileparts(mask_path);
    msk_name=['_' msk_name];
end

% Extract VOI numbers from atlas
atlas_values = unique(ATLASimg(ATLASimg>0));

% Extract Atlas values
if isempty(VOIdetails_path)
    VOIdetails_table = array2table([atlas_values,atlas_values,atlas_values]);
else
    VOIdetails_table = readtable(VOIdetails_path);
end

% Initialize table
nr_VOIS = size(VOIdetails_table,1);
nr_parameters=6;
table_out=cell(1+(nr_VOIS),1+(nr_parameters));

% Create row headers for the table
for v = 1:nr_VOIS
    if isnumeric(VOIdetails_table{v,3})
        table_out{1+v,1}=num2str(VOIdetails_table{v,3});
    else
        table_out{1+v,1}=char(VOIdetails_table{v,3}); % 1=ROIid, 3=ROIname
    end
end

% Create column headers for table
table_out{1,2} = ['mean_' img_name msk_name];
table_out{1,3} = ['median_' img_name msk_name];
table_out{1,4} = ['min_' img_name msk_name];
table_out{1,5} = ['max_' img_name msk_name];
table_out{1,6} = ['std_' img_name msk_name];
table_out{1,7} = ['n_voxels_' img_name msk_name];

% Calculate stats
colormap=zeros(size(IMG)); % initalize outcome of the segmented input image for all atlas ROIids
for v = 1:nr_VOIS
    voi=VOIdetails_table{v,1}; % grabs the ROIids
    totmask=zeros(size(IMG)); % initialize mask of the voi v
    for s=1:length(voi) % considers the case of a composite voi combining multiple ROIids, if any
        sv=voi(s); % loops through each subvoi of the (composite) voi, if any.
        if ismember(sv,atlas_values) % only consider values from the atlas image
            VOImask = zeros(size(ATLASimg));
            VOImask = VOImask + (ATLASimg == sv);
            mask = VOImask.*(MASKimg);   
            mask=mask.*sv;
            colormap=colormap+mask;  % add to colormap
            totmask=totmask+mask; % add subvoi (if any) to voi mask
        end
    end
    tmp=IMG(totmask>0);
    % Fill table
    table_out{1+v,2} = nanmean(tmp);
    table_out{1+v,3} = nanmedian(tmp);
    table_out{1+v,4} = nanmin(tmp);
    table_out{1+v,5} = nanmax(tmp);
    table_out{1+v,6} = nanstd(tmp);
    table_out{1+v,7} = sum(totmask(:)>0);
end
end



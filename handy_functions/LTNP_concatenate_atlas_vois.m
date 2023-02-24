function LTNP_concatenate_atlas_vois(atlas_path,VOItranslation_path,new_atlas_path)

% Read atlas and img
[ATLASimg,Vref] = LCN12_read_image(atlas_path);
ATLASimg = round(ATLASimg);

% Extract VOI numbers from atlas
atlas_values = unique(ATLASimg(ATLASimg>0));

% Extract Atlas values
VOItranslation_table = readtable(VOItranslation_path);

% Initialize table
nr_VOIS = size(VOItranslation_table,1);

% Calculate stats
for v = 1:nr_VOIS
    voi=VOItranslation_table{v,1}; % grabs the ROIids
    newid=VOItranslation_table{v,2};
%     for s=1:length(voi) % considers the case of a composite voi combining multiple ROIids, if any
%         sv=voi(s); % loops through each subvoi of the (composite) voi, if any.
%         if ismember(sv,atlas_values) % only consider values from the atlas image
%             ATLASimg(ATLASimg == sv)=voi;
%         end
%     end
    if ismember(voi,atlas_values) % only consider values from the atlas image
        ATLASimg(ATLASimg == voi)=newid;
    end
end

Vref.fname=new_atlas_path;
spm_write_vol(Vref,ATLASimg);

end

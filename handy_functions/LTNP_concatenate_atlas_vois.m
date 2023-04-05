function [ATLASimg,new_atlas_path]=LTNP_concatenate_atlas_vois(atlas_path,VOItranslation_path,new_atlas_path)

% Should be very carefull !
% When translating atlas values, the range of unique atlas numbers where you translate to should be
% very close to the original range of atlas numbers
% otherwise a scaling biais will be introduced when saving the image with
% spm_write_vol

% Read atlas and img
[ATLASimg,Vref] = LCN12_read_image(atlas_path);
ATLASimg = round(ATLASimg);

% Extract VOI numbers from atlas
atlas_values = unique(ATLASimg);
%atlas_values = unique(ATLASimg(ATLASimg>0));

% Extract Atlas values
VOItranslation_table = readtable(VOItranslation_path);

% Initialize table
nr_VOIS = size(VOItranslation_table,1);

% Calculate stats
for v = 1:nr_VOIS
    voi=VOItranslation_table{v,1}; % grabs the ROIids
    newid=VOItranslation_table{v,3};
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

% nr_VOIs=length(atlas_values);
% for i=1:nr_VOIs
%     if ismember(atlas_values(i),VOItranslation_table{:,1})
%         index=atlas_values(i); % given that row number corresponds to VOI number
%         ATLASimg(ismember(ATLASimg,atlas_values(i)))=VOItranslation_table{index,3};
%     end
% end


Vref.fname=new_atlas_path;
spm_write_vol(Vref,ATLASimg);

end

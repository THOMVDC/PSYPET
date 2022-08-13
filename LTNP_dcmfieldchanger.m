% change header info
% dcmfolder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/UCBJ_4s/B030';
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose';
% field_value=125860000; % has to be a number
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime';
% field_value='113945' % Means 11:39:45 and has to be a string
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife';
% field_value=1.2231e+03
% field_path='AcquisitionTime';
% field_value='164037'; % Means 16:40:37 and has to be a string
% dcmfieldchanger(dcmfolder, field_path, field_value)


function LTNP_dcmfieldchanger(dcmfolder, field_path, field_value)

T=dicomCollection(dcmfolder);

s=split(field_path,'.');

for i=1:length(T.Filenames{1})
    file=T.Filenames{1}(i);
    V = dicomread(file);
    m = dicominfo(file);
    m=setfield(m, s{:},field_value);
    dicomwrite(V, file, m, 'CreateMode', 'copy');  
end
end


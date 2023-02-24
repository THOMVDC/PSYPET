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


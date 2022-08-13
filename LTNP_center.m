function LTNP_center(input_path,output_path)

if nargin == 2
    copyfile(input_path,output_path)
    input=[output_path ',1'];
else
    input=[input_path ',1'];
end

%% Set the origin to the center of the image
% This part is written by Fumio Yamashita.
file = deblank(input); % removes trailing whitespace and null characters from a string
st.vol = spm_vol(file); % reads a volume
vs = st.vol.mat\eye(4); % A\B solves Ax=B; so it finds x so that A times x is equal to the identity matrix (eye(4)) 
vs(1:3,4) = (st.vol.dim+1)/2; % add one to your image matrix dimensions and divide by two
spm_get_space(st.vol.fname,inv(vs));

%% Set the origin to the AC by registering to mni
% Script from Carlton Chu
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810
spmDir=which('cat12');
spmDir=spmDir(1:end-7);
%tmpl=fullfile(spmDir,'templates_volumes','Template_T1_IXI555_MNI152_GS.nii'); % cat12.7
tmpl=fullfile(spmDir,'templates_1.50mm','Template_T1_IXI555_MNI152_GS.nii'); % cat12.6
vg=spm_vol(tmpl);
flags.regtype='rigid';
%p=spm_select(inf,'image');
for i=1:size(input,1)
    f=strtrim(input(i,:));
    spm_smooth(f,'temp.nii',[12 12 12]);
    vf=spm_vol('temp.nii');
    [M,~] = spm_affreg(vg,vf,flags);
    M3=M(1:3,1:3);
    [u,~,v]=svd(M3);
    M3=u*v';
    M(1:3,1:3)=M3;
    N=nifti(f);
    N.mat=M*N.mat;
    create(N);
end
end


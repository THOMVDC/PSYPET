function outputimage = LCN_3Dimage_erode(inputimage,kernel)
% inputimage = 3D binary image
% kernel = vector with 3 odd dimensions
% outputimage = 3D binary image of the same size as inputimage
% 
% Needs to be extended but for what we want this is sufficient
% author: Patrick Dupont
%         May 2013
%__________________________________________________________________________

% calculate size of the image
[dimx,dimy,dimz] = size(inputimage);
if mod(kernel(1),2) == 0
   kernel(1) = kernel(1)+1;
end
if mod(kernel(2),2) == 0
   kernel(2) = kernel(2)+1;
end
if mod(kernel(3),2) == 0
   kernel(3) = kernel(3)+1;
end

tmpx = (kernel(1)-1)/2;
tmpy = (kernel(2)-1)/2;
tmpz = (kernel(3)-1)/2;

% the bordr of the image should be removed
tmpimage = inputimage;
tmpimage(1:tmpx-1,:,:) = 0;
tmpimage(:,1:tmpy-1,:) = 0;
tmpimage(:,:,1:tmpz-1) = 0;

for i = 1:tmpx
    tmpimage(1+i:dimx,:,:) = tmpimage(i+1:dimx,:,:).*tmpimage(1:dimx-i,:,:);
end
for i = 1:tmpy
    tmpimage(:,1+i:dimy,:) = tmpimage(:,1+i:dimy,:).*tmpimage(:,1:dimy-i,:);
end
for i = 1:tmpz
    tmpimage(:,:,1+i:dimz) = tmpimage(:,:,1+i:dimz).*tmpimage(:,:,1:dimz-i);
end
outputimage = tmpimage;
end
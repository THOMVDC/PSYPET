function LTNP_apply_mask(image,mask,outname)
% image and mask should be coregistred

[IMG,Vref]=LCN12_read_image(image);
MASK=LCN12_read_image(mask);

IMGMASK=IMG.*MASK;

LCN12_write_image(IMGMASK,outname,'masked',Vref.dt(1),Vref)

end
function [answer]=LTNP_is_this_same_image(image1,image2)
[I1,Vref1]=LCN12_read_image(image1);
[I2,Vref2]=LCN12_read_image(image2);

if isequal(I1,I2) && isequal(Vref1,Vref2)
    answer=true(1);
else
    answer=false(1);
end
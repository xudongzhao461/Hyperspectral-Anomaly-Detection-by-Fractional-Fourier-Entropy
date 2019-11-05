function [FrFE,order]=FRFEorder(DataTest)
index = 1;
[rows, cols, bands]=size(DataTest);
E=zeros(11,bands);
for p = 0:0.1:1
im1 = zeros(rows, cols, bands);
for i = 1:rows
    for j = 1:cols
        im1(i,j,:) = Disfrft(squeeze(DataTest(i,j,:)),p);
    end
end
im1 = center_standard(abs(im1));

for ii = 1: bands
E(index, ii) = entropy(im1(:,:,ii));
end
index = index + 1;
end
[FrFE,order]=max(max(E,[],2));
order=order/10;
end
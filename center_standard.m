function [ output ] = center_standard( input )
%CENTER_STANDARD 此处显示有关此函数的摘要
%   此处显示详细说明
temp = input(:);
sigma = cov(temp);
sigma = sqrt(sigma);
in_mean = mean(temp);
output = (input - in_mean)/sigma;
end


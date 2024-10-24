function [out] = rowNormAccrossM(A)
%ROWNORMACCROSSM Summary of this function goes here
%   Detailed explanation goes here
out=sqrt(A(:,:,1).^2 + A(:,:,2).^2 + A(:,:,3).^2);
end


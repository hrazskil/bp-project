function [directivity] = directivity(radiativity,powerFull)
%DIRECTIVITY Summary of this function goes here
%   Detailed explanation goes here
directivity = (4*pi*radiativity)./powerFull;
end


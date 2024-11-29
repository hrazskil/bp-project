function [x,y,z] = sph2cartKH(r,theta,phi)
%CART2SPHKH Summary of this function goes here
%Detailed explanation goes here
x       = r.*cos(phi).*sin(theta);
y       = r.*sin(phi).*sin(theta);
z       = r.*cos(theta);
end


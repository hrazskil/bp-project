function [r,theta,phi] = cart2sphKH(x,y,z)
%CART2SPHKH Summary of this function goes here
%Detailed explanation goes here
r               = sqrt(x.^2+y.^2+z.^2);
theta           = angle(z + 1i*sqrt(x.^2+y.^2));
theta(isnan(theta)) = 0;
phi             = angle(x + 1i*y);
phi(isnan(phi)) = 0;
end


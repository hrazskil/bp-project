function [Fr,Fth,Fph] = cart2sphVec(x,y,z, Fx, Fy, Fz)
%CART2SPH0 Summary of this function goes here
%   Detailed explanation goes here

% repeated calculations
r       = sqrt(x.^2 + y.^2 + z.^2);
xy      = sqrt(x.^2 + y.^2);

% Calculate
Fr = (x.*Fx + y.*Fy + z.*Fz)./r;
Fth = (Fx.*z.*x + Fy.*z.*y -Fz.*xy.^2)./(r.*xy);
Fph = (-Fx.*y + Fy.*x)./xy;
end
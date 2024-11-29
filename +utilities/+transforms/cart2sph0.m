function [r0,theta0,phi0] = cart2sph0(x,y,z)
%CART2SPH0 Summary of this function goes here
%   Detailed explanation goes here

% repeated calculations
r       = sqrt(x.^2 + y.^2 + z.^2);
xy      = sqrt(x.^2 + y.^2);

% Initialize matrices
r0      = zeros(size(x, 1), 3);
theta0  = zeros(size(x, 1), 3);
phi0    = zeros(size(x, 1), 3);

% Calculate
r0(:,1) = (x)./r;
r0(:,2) = (y)./r;
r0(:,3) = (z)./r;
r0(isnan(r0)) = 0;

theta0(:,1) = (z.*x)./(r.*xy);
theta0(:,2) = (z.*y)./(r.*xy);
theta0(:,3) = -(x.^2+y.^2)./(r.*xy);
theta0(isnan(theta0)) = 0;

phi0(:,1) = -y./xy;
phi0(:,2) = x./xy;
phi0(:,3) = 0;
phi0(isnan(phi0)) = 0;
end
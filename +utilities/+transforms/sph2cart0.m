function [x0,y0,z0] = sph2cart0(theta,phi)
%SPH2CART0 Summary of this function goes here
%   Detailed explanation goes here
sthe    = sin(theta);
cthe    = cos(theta);
sphi    = sin(phi);
cphi    = cos(phi);

% Initialize matrices
stheta = size(theta, 1);
x0       = zeros(stheta, 3);
y0       = zeros(stheta, 3);
z0       = zeros(stheta, 3);

% Calculate
x0(:,1) = sthe.*cphi;
x0(:,2) = cthe.*cphi;
x0(:,3) = -sphi;


y0(:,1) = sthe.*sphi;
y0(:,2) = cthe.*sphi;
y0(:,3) = cphi;


z0(:,1) = cthe;
z0(:,2) = -sthe;
z0(:,3) = 0;

end


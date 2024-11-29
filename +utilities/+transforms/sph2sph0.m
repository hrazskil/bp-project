function [r0,theta0,phi0] = sph2sph0(theta,phi)
%SPH2SPH0 Summary of this function goes here
%   Detailed explanation goes here

% Initialize matrices
stheta = size(theta, 1);
r0      = zeros(stheta, 3);
theta0  = zeros(stheta, 3);
phi0    = zeros(stheta, 3);

sthe    = sin(theta);
cthe    = cos(theta);

sphi    = sin(phi);
cphi    = cos(phi);
% Calculate
r0(:,1) = sthe.*cphi;
r0(:,2) = sthe.*sphi;
r0(:,3) = cthe;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

theta0(:,1) = cthe.*cphi;
theta0(:,2) = cthe.*sphi;
theta0(:,3) = -sthe;

theta0(theta0 == 1) = 0;

phi0(:,1) = -sphi;
phi0(:,2) = cphi;
phi0(:,3) = 0;
end


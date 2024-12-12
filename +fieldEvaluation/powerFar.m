function [pow] = powerFar(eF,rObserved)
%POWERFAR Summary of this function goes here
%   Detailed explanation goes here
const   = utilities.constants.giveConstants;
[r0]   = utilities.transforms.cart2sph0( ...
    rObserved(:,1),rObserved(:,2),rObserved(:,3))
pow     = (1/(2*const.Z0))*(abs(eF).^2)*r0;
end


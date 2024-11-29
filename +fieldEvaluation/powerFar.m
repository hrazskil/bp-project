function [pow] = powerFar(eF)
%POWERFAR Summary of this function goes here
%   Detailed explanation goes here
const = utilities.constants.giveConstants;
pow = (1/(2*const.Z0))*(abs(eF).^2);
end


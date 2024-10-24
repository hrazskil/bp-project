function [pow] = powerFar(eleField)
%POWERFAR Summary of this function goes here
%   Detailed explanation goes here
const.utilities.constants.giveConstants;
pow = 1/(2*const.Z0)*(abs(eleField).^2);
end


function [pow] = powerFar(eF, rObserved)
%POWERFAR Computes the power in the far-field from the electric field and observed points.
%   Inputs:
%       eF - Nx3 matrix of electric field vectors at observed points.
%       rObserved - Nx3 matrix of observed positions (x, y, z).
%   Output:
%       pow - Power density at observed points.

const = utilities.constants.giveConstants; % Retrieve constants
r0 = utilities.rowNorm(rObserved); % Compute the radial distance to each observation point

% Calculate the power density using the magnitude squared of E
pow = (1 / (2 * const.Z0)) * (abs(eF).^2) .* r0.^2; % r0^2 accounts for the distance-squared
end
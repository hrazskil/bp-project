function [powerDensity] = powerDensityFar(eF)
% powerDensityFar Computes far-field power density from electric field vectors.
%
% INPUTS:
%   eF (nObs × 3) - Complex electric field vectors at observation points (V/m).
%
% OUTPUT:
%   powerDensity (nObs × 1) - Power density at each point (W/m²).
%
% FORMULATION:
%   S = |E|² / (2 * Z0)
%   where Z0 is the free-space impedance, and E is the far-field electric vector.
%
% ------------------------------------------------------------------------

% Step 1: Retrieve physical constant
const = utilities.constants.giveConstants;
Z0 = const.Z0;  % Free-space impedance (Ω)

% Step 2: Compute power density from electric field magnitude
powerDensity = (1 / (2 * Z0)) * sum(abs(eF).^2, 2); % nObs × 1

end

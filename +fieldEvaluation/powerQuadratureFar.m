function [integral] = powerQuadratureFar(Nleb, dip, f)
% powerQuadratureFar Estimates radiated power using far-field approximation.
%
% INPUTS:
%   Nleb (scalar) - Number of Lebedev quadrature points.
%   dip  (struct) - Dipole configuration (positions, amplitudes, directions).
%   f    (scalar) - Frequency in Hz.
%
% OUTPUT:
%   integral (scalar) - Total radiated power estimate (W).
%
% FORMULATION:
%   Uses far-field |F|² scaled by 1 / (2 * Z0) over the unit sphere.
%
% ------------------------------------------------------------------------

% Step 1: Get quadrature points and weights
[points, weights, ~] = utilities.getLebedevSphere(Nleb);

% Step 2: Physical constants
construct = utilities.constants.giveConstants;

% Step 3: Evaluate far-field electric field
fF = fieldEvaluation.farFieldM2(points, dip, f);

% Step 4: Compute integrated power
integral = sum(sum(fF .* conj(fF), 2) .* weights) / (2 * construct.Z0);

end

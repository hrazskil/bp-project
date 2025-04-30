function [propagationFactor] = computePropagationFactor(rObserved, k)
% computePropagationFactor Computes the far-field propagation factor.
%
% INPUTS:
%   rObserved (nObs × 3) - Cartesian coordinates of the observation points.
%   k         (scalar)   - Wavenumber (2πf / c).
%
% OUTPUT:
%   propagationFactor (nObs × 1) - Complex amplitude scaling factor accounting for
%                                  phase delay and 1/r decay.
%
% FORMULATION:
%   propagationFactor = exp(-1i * k * r) / r
%   where r is the Euclidean distance to each observation point.
%
% ------------------------------------------------------------------------

% Step 1: Compute distance from origin to each observation point
r = utilities.rowNorm(rObserved); % Euclidean norms (nObs × 1)

% Step 2: Compute complex propagation factor
propagationFactor = exp(-1i * k * r) ./ r;

end
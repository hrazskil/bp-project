function [kVec] = computeWaveVector(rObserved, dipPos, freq)
% computeWaveVector Computes wave vectors from dipole positions to observation points.
%
% INPUTS:
%   rObserved (nObs × 3) - Cartesian observation points.
%   dipPos    (nObs × 3) - Corresponding dipole positions.
%   freq      (scalar)   - Frequency in Hz.
%
% OUTPUT:
%   kVec (nObs × 3) - Wave vectors (direction × magnitude).
%
% ------------------------------------------------------------------------

% Step 1: Get constants and compute wavenumber
construct = utilities.constants.giveConstants;
omega = 2 * pi * freq;
k_mag = omega / construct.c0;

% Step 2: Relative vectors and normalization
rVec = rObserved - dipPos; % Direction vectors
k_hat = rVec ./ utilities.rowNorm(rVec); % Unit vectors

% Step 3: Apply magnitude
kVec = k_mag * k_hat;

end

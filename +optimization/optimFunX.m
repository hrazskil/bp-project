function out = optimFunX(amp, dip, inputData, numDipoles)
% optimFunX Objective function wrapper for PSO optimization of dipole amplitudes.
%
%   This function reconstructs the complex dipole amplitudes from the
%   input real-valued vector amp and evaluates the normalized error
%   between the simulated and reference far-field patterns.
%
% INPUTS:
%   amp        (2 * numDipoles × 1) - Real-valued vector of concatenated 
%                                     real and imaginary parts of dipole amplitudes.
%   dip        (struct)             - Dipole structure containing:
%       dip.pos        (numDipoles × 3) - Dipole positions.
%       dip.dir        (numDipoles × 3) - Dipole orientations (unit vectors).
%   inputData  (struct)             - Structure with far-field reference data and setup:
%       inputData.rObserved         - Observation points.
%       inputData.referenceFields   - Reference far-field values.
%       inputData.weights           - Optional weighting for error computation.
%   numDipoles (scalar)             - Total number of dipoles.
%
% OUTPUT:
%   out (scalar) - Objective function value (normalized error).
%
% --------------------------------------------------------------------------------------------

%% Step 1: Reconstruct Complex Dipole Amplitudes
% Extract real and imaginary parts from amp vector
realPart = amp(1:numDipoles);
imagPart = amp(numDipoles+1:2*numDipoles);

% Combine into complex amplitudes
dip.complAmpl = realPart.' + 1i * imagPart.';

%% Step 2: Evaluate Objective Function
% Compute the normalized error between simulated and reference far-fields
out = optimization.normObjectiveFunction_rad(dip, inputData);

end

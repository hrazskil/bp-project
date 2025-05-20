function [fF] = farFieldM2(rObserved, dip, f)
% farFieldM2 Computes the far-field electric radiation pattern from a set of dipoles.
%
% INPUTS:
%   rObserved (nObs × 3)  - Cartesian coordinates of the observation points.
%   dip        (struct)   - Structure containing dipole properties:
%       dip.pos (nDip × 3)        - Cartesian positions of the dipoles.
%       dip.complAmpl (nDip × 1)  - Complex amplitudes of the dipoles.
%       dip.dir (nDip × 3)        - Dipole moment directions (unit vectors).
%   f          (scalar)   - Frequency of the electromagnetic wave in Hz.
%
% OUTPUT:
%   fF (nObs × 3) - Far-field electric field vectors at each observation point.
%
% FORMULATION:
% The far-field approximation is used to compute the radiated electric field from each dipole.
% The total field is obtained by superposing the contributions from all dipoles.
%
% --------------------------------------------------------------------------------------------

%% Step 1: Initialization
nDip = size(dip.pos, 1);    % Number of dipoles
nObs = size(rObserved, 1);  % Number of observation points

%% Step 2: Define Physical Constants
construct = utilities.constants.giveConstants; % Retrieve fundamental constants
omega = 2 * pi * f;         % Angular frequency
k = omega / construct.c0;   % Wavenumber (k = omega / c)

%% Step 3: Compute Unit Vectors in Observation Directions
% Normalize the observation vectors to get direction unit vectors (r̂)
rHat = rObserved ./ repmat(utilities.rowNorm(rObserved), [1, 3]); 
rHat = repmat(rHat, [nDip, 1]);  % Expand to match dipole count

%% Step 4: Compute Dipole Moments
% Expand dipole amplitudes and directions to match each observation point
Mp = repelem(dip.complAmpl, nObs, 1) .* repelem(dip.dir, nObs, 1);

%% Step 5: Compute Phase Term (Far-field Approximation)
% Phase factor: exp(i * k * r̂ ⋅ r_dip)
phase = exp(1i * k * sum(rHat .* repelem(dip.pos, nObs, 1), 2));

%% Step 6: Compute Electric Field Using Dipole Far-Field Expression
% Vector cross product formulation for far-field electric field:
% E_far = (Z0 * c * k^2 / 4π) * exp(i * k * r̂ ⋅ r_dip) * [-r̂ × (r̂ × Mp)]
fF = construct.Z0 * construct.c0 * k^2 / (4 * pi) * ...
     (phase .* utilities.rowCross(-rHat, utilities.rowCross(rHat, Mp)));

%% Step 7: Sum Contributions from All Dipoles
% Reshape and sum along the dipole dimension to get total field at each observation point
fF = reshape(sum(reshape(fF.', [], nDip), 2), 3, [])';

end

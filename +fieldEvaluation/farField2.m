function [fF] = farField2(rObserved, dip, f)
% FARFIELD Computes the far-field radiation pattern from an array of dipoles.
%
% This function calculates the far-field approximation of the electric field 
% radiated by a set of dipoles at given observation points.
%
% INPUTS:
%   rObserved (nObs × 3) - Cartesian coordinates of observation points.
%   dip       (struct)   - Structure containing dipole properties:
%       dip.pos (nDip × 3)        - Positions of the dipoles.
%       dip.complAmpl (nDip × 1)  - Complex amplitudes of the dipoles.
%       dip.dir (nDip × 3)        - Dipole moment directions (unit vectors).
%   f         (scalar)   - Frequency of the electromagnetic wave in Hz.
%
% OUTPUT:
%   fF (nObs × 3) - Far-field electric field vectors at each observation point.
%
% ------------------------------------------------------------------------

%% **Step 1: Initialization**
nDip = size(dip.pos, 1);    % Number of dipoles
nObs = size(rObserved, 1);  % Number of observation points

%% **Step 2: Define Physical Constants**
construct = utilities.constants.giveConstants; % Retrieve fundamental constants
omega = 2 * pi * f;         % Angular frequency
k = omega / construct.c0;   % Wavenumber (k = 2π / λ)

%% **Step 3: Normalize Observation Point Directions**
% Normalize each observed position vector to get unit direction vectors
rNorms = utilities.rowNorm(rObserved);   % Compute norms of observation points
r0 = rObserved ./ repmat(rNorms, [1, 3]); % Normalize to unit vectors

% Expand r0 to match the number of dipoles (each dipole contributes to all points)
r0 = repmat(r0, [nDip, 1]); 

%% **Step 4: Compute Dipole Moments**
% Multiply each dipole's amplitude and direction, then repeat for all observation points
Mp = repelem(dip.complAmpl, nObs, 1) .* repelem(dip.dir, nObs, 1);

%% **Step 5: Compute Far-Field Contribution**
% The phase factor for each dipole at each observation point
phaseFactor = exp(1i * k * sum(r0 .* repelem(dip.pos, nObs, 1), 2));

% Compute the far-field electric field components
fF = construct.Z0 * construct.c0 * k^2 .* phaseFactor ./ (4 * pi) .* ...
    (utilities.rowCross(-r0, utilities.rowCross(r0, Mp)));

%% **Step 6: Sum Contributions from All Dipoles**
% Reshape and sum along the dipole dimension to get the total field at each observation point
fF = reshape(sum(reshape(fF.', [], nDip), 2), 3, [])'; 

end

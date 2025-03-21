function [eM] = magFieldM3(rObserved, dip, f)
% magFieldM3 Computes the magnetic field at given observation points due to a set of dipoles.
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
%   eM (nObs × 3) - Magnetic field vectors at each observation point.
%
% FORMULATION:
% The magnetic field is computed using the dipole radiation formula.
% The function sums up contributions from all dipoles at each observation point.
%
% ------------------------------------------------------------------------

%% Step 1: Initialization
nDip = size(dip.pos, 1);    % Number of dipoles
nObs = size(rObserved, 1);  % Number of observation points

%% Step 2: Define Physical Constants
construct = utilities.constants.giveConstants; % Retrieve fundamental constants
omega = 2 * pi * f;         % Angular frequency
k = omega / construct.c0;   % Wavenumber (k = 2π / wavelength)

%% Step 3: Compute Distance Vectors Between Dipoles and Observation Points
% Compute relative position vectors from dipoles to observation points
MRVec = repmat(rObserved, nDip, 1) - repelem(dip.pos, nObs, 1); 

% Compute the Euclidean distance between dipoles and observation points
MR = utilities.rowNorm(MRVec); % MR is an (nDip * nObs) × 1 vector

% Normalize distance vectors (unit vectors pointing from dipoles to observation points)
MdirR = MRVec ./ MR; 

%% Step 4: Compute Dipole Moments in Space
% Expand dipole complex amplitudes and orientations for each observation point
Mp = repelem(dip.complAmpl, nObs, 1) .* repelem(dip.dir, nObs, 1);

%% Step 5: Compute Magnetic Field Components
% Compute the phase term (1i * k * R)
IMR = 1i * k * MR;

% Compute the magnetic field using the dipole radiation formula
eM = construct.c0 * k^2 * exp(-IMR) ./ (4 * pi * MR) .* ...
    (utilities.rowCross(MdirR, Mp) .* (1 ./ IMR + 1));

%% Step 6: Sum Contributions from All Dipoles
% Reshape and sum along the dipole dimension to get total field at each observation point
eM = reshape(sum(reshape(eM', [], nDip), 2), 3, [])';

end
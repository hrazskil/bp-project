function dip = halfwaveDipoleArray(dipoleLength, numDipoles, rotationAngleRad, positionOffset, amplitudes)
%HALFWAVEDIPOLEARRAY Constructs a discretized half-wave dipole using elementary dipoles.
%
%   Inputs:
%       dipoleLength      - Total length of the physical dipole [meters]
%       numDipoles        - Number of elementary dipoles (>=2)
%       rotationAngleRad  - Rotation angle around x-axis [radians]
%       positionOffset    - 1×3 translation vector [meters]
%       amplitudes        - Optional complex amplitudes (numDipoles×1)
%
%   Output:
%       dip - Struct with fields:
%           .pos        - (numDipoles)×3 dipole center positions
%           .dir        - (numDipoles)×3 unit direction vectors
%           .complAmpl  - (numDipoles)×1 complex amplitudes

    if numDipoles < 2
        error('numDipoles must be at least 2 to span the dipole length.');
    end

    % Generate z-coordinates of dipole centers (centered around zero)
    zCenters = linspace(-dipoleLength/2, dipoleLength/2, numDipoles).';

    % Rotation matrix around x-axis
    Rx = [1, 0, 0;
          0, cos(rotationAngleRad), -sin(rotationAngleRad);
          0, sin(rotationAngleRad),  cos(rotationAngleRad)];

    % Initial positions along z-axis
    centers = [zeros(numDipoles,1), zeros(numDipoles,1), zCenters];

    % Apply rotation and translation
    centers = (Rx * centers.').';
    centers = centers + positionOffset;

    % Default direction (along z) and rotate
    dirs = repmat([0, 0, 1], numDipoles, 1);
    dirs = (Rx * dirs.').';

    % Default amplitudes if none given
    if nargin < 5 || isempty(amplitudes)
        amplitudes = ones(numDipoles, 1);
    end

    if length(amplitudes) ~= numDipoles
        error('Length of amplitudes must match numDipoles.');
    end

    % Assemble output
    dip.pos = centers;
    dip.dir = dirs;
    dip.complAmpl = amplitudes;
end

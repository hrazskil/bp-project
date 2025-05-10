% Parameters
numPairs = 4;
dipoleLength = 0.05;                 % Total dipole length [m]
numDipoles = 50;
offsetStep = 0.05;                   % Offset between dipole pairs [m]
rotationY = pi/4;                    % Rotation for y-direction dipoles
rotationX = 3*pi/4;                  % Rotation for x-direction dipoles


% Optional: define excitation amplitudes
amplitudes = ones(numDipoles, 1);  % Uniform excitation

% Initialize empty structure for full array
dip.pos = [];
dip.dir = [];
dip.complAmpl = [];

% Generate crossed dipole pairs
for i = 1:numPairs
    % Vertical z-offset for stacking
    zOffset = (i - 1) * offsetStep;
    posOffset = [0, 0, zOffset];

    % Generate dipole rotated for Y orientation
    dipY = geometry.halfwaveDipoleArray(dipoleLength, numDipoles, rotationY, posOffset, amplitudes);

    % Generate dipole rotated for X orientation
    dipX = geometry.halfwaveDipoleArray(dipoleLength, numDipoles, rotationX, posOffset, amplitudes);

    % Append to full structure
    dip.pos        = [dip.pos; dipY.pos; dipX.pos];
    dip.dir        = [dip.dir; dipY.dir; dipX.dir];
    dip.complAmpl  = [dip.complAmpl; dipY.complAmpl; dipX.complAmpl];
end

% === Optional Visualization ===
figure;
quiver3(dip.pos(:,1), dip.pos(:,2), dip.pos(:,3), ...
        dip.dir(:,1), dip.dir(:,2), dip.dir(:,3), ...
        0.01, 'LineWidth', 1.2);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Crossed Half-Wave Dipole Pairs (Stacked)');
grid on;

% === Save dipole structure ===
save('dipoleStructure.mat', 'dip');
disp('Dipole structure saved to halfwaveDipole.mat.');
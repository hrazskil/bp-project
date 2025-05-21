% === Parameters for Crossed Dipole "Pluses" ===
numPluses = 4;                    % Total number of pluses
dipoleLength = 1;              % Length of each dipole [m]
numDipolesPerArm = 20;            % Dipoles per arm (X and Z)
offsetStep = 1.5;                % Spacing between pluses along Z
rotationX = pi/2;                 % Dipoles along Z (vertical arm)
rotationY = 0;                    % Dipoles along X (horizontal arm)

amplitudes = ones(numDipolesPerArm, 1);  % Uniform excitation

% === Initialize dipole structure ===
dip.pos = [];
dip.dir = [];
dip.complAmpl = [];

% === Generate Crossed Pluses ===
for i = 1:numPluses
    % Centered around Z = 0
    zOffset = (i - (numPluses + 1)/2) * offsetStep;
    origin = [0, 0, zOffset];  % All at x=0, y=0

    % Horizontal dipole line (X-direction)
    dipX = geometry.halfwaveDipoleArray(dipoleLength, numDipolesPerArm, rotationY, origin, amplitudes);

    % Vertical dipole line (Z-direction)
    dipZ = geometry.halfwaveDipoleArray(dipoleLength, numDipolesPerArm, rotationX, origin, amplitudes);

    % Append both arms to full dipole array
    dip.pos        = [dip.pos; dipX.pos; dipZ.pos];
    dip.dir        = [dip.dir; dipX.dir; dipZ.dir];
    dip.complAmpl  = [dip.complAmpl; dipX.complAmpl; dipZ.complAmpl];
end

% === Grid Parameters (YZ-plane) ===
gridSizeY = 5;                   % Number of points along Y
gridSizeZ = 5;                   % Number of points along Z
gridWidthY = 1;                % Total physical width along Y [m]
gridHeightZ = 1;               % Total physical height along Z [m]
gridOffsetX = -0.38;             % X-offset behind each plus
rotationGrid = 0;                % Dipoles oriented along X

% Compute spacing from desired total width and height
gridSpacingY = gridWidthY / (gridSizeY - 1);
gridSpacingZ = gridHeightZ / (gridSizeZ - 1);

% Amplitudes for each grid dipole
gridDipolesPerPlane = gridSizeY * gridSizeZ;
amplitudesGrid = ones(gridDipolesPerPlane, 1);

% Optional: compute height of plus to align grid
plusHeight = (numDipolesPerArm - 1) * dipoleLength;

for i = 1:numPluses
    % Z position of the plus
    zOffset = (i - (numPluses + 1)/2) * offsetStep;

    % Centered in Y and Z to match the plus geometry
    [Yg, Zg] = meshgrid( ...
        linspace(-gridWidthY/2, gridWidthY/2, gridSizeY), ...
        linspace(zOffset - gridHeightZ/2, zOffset + gridHeightZ/2, gridSizeZ) ...
    );
    Xg = ones(size(Yg)) * gridOffsetX;
    gridPos = [Xg(:), Yg(:), Zg(:)];

    % Generate dipoles
    dipGrid = geometry.halfwaveDipoleArray(dipoleLength, gridDipolesPerPlane, ...
                 rotationGrid, [0, 0, 0], amplitudesGrid);
    dipGrid.pos = gridPos;

    % Append to dipole structure
    dip.pos        = [dip.pos; dipGrid.pos];
    dip.dir        = [dip.dir; dipGrid.dir];
    dip.complAmpl  = [dip.complAmpl; dipGrid.complAmpl];
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
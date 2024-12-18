% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%%
% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\test\halfwaveDipole.mat');

% Define the number of theta and phi points
nTh = 50; 
nPh = 99;

% Create linearly spaced vectors for theta and phi
theta = linspace(0, pi, nTh); 
phi = linspace(0, 2*pi, nPh);

% Create a meshgrid for theta and phi
[Theta, Phi] = meshgrid(theta, phi);

% Reshape Theta and Phi into column vectors
Theta = reshape(Theta, [nTh*nPh, 1]);
Phi = reshape(Phi, [nTh*nPh, 1]);

% Convert spherical coordinates to Cartesian coordinates
[x, y, z] = utilities.transforms.sph2cartKH(ones(nTh*nPh, 1), Theta, Phi);
rObserved = [x, y, z]; % Combine Cartesian coordinates into a matrix

% Evaluate the far field based on the observed coordinates
[fF] = fieldEvaluation.farField(rObserved, dip, f0List);

% Transform Cartesian field components to spherical coordinates
[Fr, Fth, Fph] = utilities.transforms.cart2sphVec(...
    x, y, z, fF(:, 1), fF(:, 2), fF(:, 3) ...
    );

% Reshape Theta, Phi, Fth, and Fph for contour plotting
Theta = reshape(Theta, [nPh, nTh]);
Phi = reshape(Phi, [nPh, nTh]);
Fth = reshape(Fth, [nPh, nTh]);
Fph = reshape(Fph, [nPh, nTh]);

% Visualization of Fth and Fph using contour plots with color bars

% Plot for the real part of Fth
figure
contourf(Theta/pi, Phi/(2*pi), real(Fth));
title('Re[Fth] ver2');
colorbar; % Add color bar for reference
grid on;

% Plot for the imaginary part of Fth
figure
contourf(Theta/pi, Phi/(2*pi), imag(Fth));
title('Im[Fth] ver2');
colorbar; 
grid on;

% Plot for the real part of Fph
figure
contourf(Theta/pi, Phi/(2*pi), real(Fph));
title('Re[Fph] ver2');
colorbar; 
grid on;

% Plot for the imaginary part of Fph
figure
contourf(Theta/pi, Phi/(2*pi), imag(Fph));
title('Im[Fph] ver2');
colorbar;
grid on;
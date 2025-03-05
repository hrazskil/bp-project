% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%%
% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\halfwaveDipole.mat');

% Define the number of theta and phi points
nTh = 120; 
nPh = 240;

% Create linearly spaced vectors
theta = linspace(0, pi, nTh); 
phi = linspace(0, 2*pi, nPh);

% Create a meshgrid
[Phi, Theta] = meshgrid(phi,theta);

% Reshape into column vectors
Theta = Theta(:); % Flatten into a single column vector
Phi = Phi(:); % <=> Phi = reshape(Phi, [nTh*nPh, 1]);


% Convert spherical to Cartesian
[x, y, z] = utilities.transforms.sph2cartCoor(ones(nTh*nPh, 1), Theta, Phi);
rObserved = [x, y, z]; % Combine Cartesian coordinates into a matrix

% Evaluate the far field
[fF] = fieldEvaluation.farField(rObserved, dip, f0List);

% Transform Cartesian field components to spherical
[Fr, Fth, Fph] = utilities.transforms.cart2sphVec(...
    x, y, z, fF(:, 1), fF(:, 2), fF(:, 3) ...
                                                 );


% Reshape Theta, Phi, Fth, and Fph for contour plotting
Theta = reshape(Theta, [nTh, nPh]);
Phi = reshape(Phi, [nTh, nPh]);
Fth = reshape(Fth, [nTh, nPh]);
Fph = reshape(Fph, [nTh, nPh]);

% Visualization of Fth and Fph using contour plots with color bars

% Plot for the real part of Fth
figure
contourf(Theta/pi, Phi/(pi), real(Fth));
title('Re[Eth]');
colorbar; % Add color bar for reference


% Plot for the imaginary part of Fth
figure
contourf(Theta/(pi), Phi/pi, imag(Fth));
title('Im[Eth]');
colorbar; 


% Plot for the real part of Fph
figure
contourf(Theta/pi, Phi/(pi), real(Fph));
title('Re[Eph]');
colorbar; 


% Plot for the imaginary part of Fph
figure
contourf(Theta/pi, Phi/(pi), imag(Fph));
title('Im[Eph]');
colorbar;

%%

% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\halfwaveDipoleFarFields.mat')

theta = linspace(0, pi, 60);
phi = linspace(0, 2*pi, 120);
[PH,TH] = meshgrid(phi,theta);
figure
contourf(TH/pi,PH/pi,real(farfield.FTheta))
xlabel('\theta/ \pi')
ylabel('\phi/ \pi')
title('Re[Fth]')

figure
pcolor(abs(farfield.FTheta -Fth))
abs(farfield.FTheta(30,60))
abs(Fth(30,60))
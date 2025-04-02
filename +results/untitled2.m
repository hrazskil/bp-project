%% --- 1_a Initialization ---
clc; clear; close all;

%% directivity vectors extraction
% Load dipole variables from file 
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\test_structure_1\halfwaveDipole.mat')

%% --- 1_c Compute Far-Field Parameters (Uniform Grid) ---
% Define the number of points for theta and phi
Ntheta = 360;  % Number of points for theta (0 to pi)
Nphi = 360;    % Number of points for phi (0 to 2pi)

% Define the angle values for theta and phi
theta_vals = linspace(0, 2*pi, Ntheta).';  % Theta values from 0 to pi
phi_vals = linspace(0, 2*pi, Nphi).';    % Phi values from 0 to 2pi

% Create grid of theta and phi values for the planes
[theta_grid1, phi_grid1] = meshgrid(theta_vals, 0);  % For Plane 1 (theta varies, phi = 0)
[theta_grid2, phi_grid2] = meshgrid(pi/2, phi_vals);  % For Plane 2 (phi varies, theta = pi/2)

% Convert spherical to Cartesian coordinates for Plane 1 (theta varies, phi = 0)
r_1 = ones(Ntheta, 1);  % Fixed distance for simplicity
[x1, y1, z1] = utilities.transforms.sph2cartCoor(r_1, theta_grid1.', phi_grid1.');  % Plane 1 points

% Convert spherical to Cartesian coordinates for Plane 2 (phi varies, theta = pi/2)
r_2 = ones(Nphi, 1);  % Fixed distance for simplicity
[x2, y2, z2] = utilities.transforms.sph2cartCoor(r_2, theta_grid2, phi_grid2);  % Plane 2 points

% Now, x1, y1, z1 contain the cartesian coordinates for Plane 1
% x2, y2, z2 contains the cartesian coordinates for Plane 2

% combine these points into a single matrix
vertical_points = [x1(:), y1(:), z1(:)];  % Flattening to get a list of points
horizontal_points = [x2(:), y2(:), z2(:)];  % Flattening to get a list of points


% Compute the far-field at the observation points
farField_plane1 = fieldEvaluation.farField(vertical_points, dip, f0List);  % Example function call
farField_plane2 = fieldEvaluation.farField(horizontal_points, dip, f0List);  % Example function call

rad_Vertical = sum(abs(farField_plane1).^2, 2);
rad_Horizontal = sum(abs(farField_plane2).^2, 2);

data = [theta_vals, rad_Vertical];
save('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\test_structure_1\dataVertical.tsv','data','-ascii')
data = [phi_vals, rad_Horizontal];
save('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\test_structure_1\dataHorizontal.tsv','data','-ascii')

%% --- 1_b Initialization 2 ---
clc; clear; close all;

%% directivity vectors extraction
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\test_structure_1\dataVertical.tsv')

load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\test_structure_1\dataHorizontal.tsv')

load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\test_structure_1\halfwaveDipole.mat')

% Create grid of theta and phi values for the planes
[theta_grid1, phi_grid1] = meshgrid(dataHorizontal(:,1), 0);  % For Plane 1 (theta varies, phi = 0)
[theta_grid2, phi_grid2] = meshgrid(pi/2, dataVertical(:,1));  % For Plane 2 (phi varies, theta = pi/2)

Ntheta = size(dataHorizontal(:,2),1);  % Number of points for theta (0 to pi)
Nphi = size(dataVertical(:,2),1);    % Number of points for phi (0 to 2pi)

% Convert spherical to Cartesian coordinates for Plane 1 (theta varies, phi = 0)
r_1 = ones(Ntheta, 1);  % Fixed distance for simplicity
[x1, y1, z1] = utilities.transforms.sph2cartCoor(r_1, theta_grid1.', phi_grid1.');  % Plane 1 points

% Convert spherical to Cartesian coordinates for Plane 2 (phi varies, theta = pi/2)
r_2 = ones(Nphi, 1);  % Fixed distance for simplicity
[x2, y2, z2] = utilities.transforms.sph2cartCoor(r_2, theta_grid2, phi_grid2);  % Plane 2 points

% combine these points into a single matrix
vertical_points = [x1(:), y1(:), z1(:)];  % Flattening to get a list of points
horizontal_points = [x2(:), y2(:), z2(:)];  % Flattening to get a list of points

%% defining input data structure
inputData.horizontal.rad        = dataHorizontal(:,2);
inputData.horizontal.points     = horizontal_points;
inputData.horizontal.weights    = ones(Ntheta, 1);

% sizeStep = dataHorizontal(2)-dataHorizontal(1);
% totalPower = sum((inputData.horizontal.rad.*inputData.horizontal.weights).*sizeStep);
totalPower_hor = trapz(inputData.horizontal.rad.*inputData.horizontal.weights);

inputData.horizontal.totalPower = totalPower_hor;

    
inputData.vertical.rad          = dataVertical(:,2);
inputData.vertical.points       = vertical_points;
inputData.vertical.weights      = ones(Ntheta,1);

% sizeStep = dataVertical(2)-dataVertical(1);
% totalPower = sum((inputData.vertical.rad.*inputData.vertical.weights).*sizeStep);
totalPower_ver = trapz(inputData.vertical.rad.*inputData.vertical.weights);

inputData.vertical.totalPower = totalPower_ver;

inputData.freq                  = f0List;        

%% --- 1_e Validate Objective Function Before Optimization ---
ampGuess = dip.complAmpl; % dipolePerturbedRef.complAmpl;
dip.complAmpl = ampGuess; % Ensure ampGuess in dip
error = optimization.normObjectiveFunction_rad(dip, inputData);
disp(['Initial Test Error: ', num2str(error)]);


%% --- Normalize Reference Dipole Amplitudes ---
numDipoles = numel(dip.complAmpl);
dipoleAmpReal = real(dip.complAmpl);
dipoleAmpImag = imag(dip.complAmpl);

% Avoid division by zero
epsilon = 1e-12;
maxAmpReal = max(abs(dipoleAmpReal), epsilon);
maxAmpImag = max(abs(dipoleAmpImag), epsilon);

normalizedAmpReal = dipoleAmpReal / maxAmpReal;
normalizedAmpImag = dipoleAmpImag / maxAmpImag;

% Save reference dipole amplitudes after normalization
dipoleRef = dip;
dipoleRef.complAmpl = normalizedAmpReal + 1i * normalizedAmpImag;

%% --- 1_c Perturbation of Initial Amplitudes ---
realPerturbationFactor = 1 + 0.01 * randn(numDipoles, 1);
imagPerturbationFactor = 1 + 0.01 * randn(numDipoles, 1);

% Apply perturbations separately to real and imaginary parts
perturbedAmp = real(dip.complAmpl) .* realPerturbationFactor + ...
               1i * imag(dip.complAmpl) .* imagPerturbationFactor;

% Normalize perturbed dipole amplitudes
perturbedAmpReal = real(perturbedAmp);
perturbedAmpImag = imag(perturbedAmp);

normalizedPerturbedAmpReal = perturbedAmpReal / max(abs(perturbedAmpReal), epsilon);
normalizedPerturbedAmpImag = perturbedAmpImag / max(abs(perturbedAmpImag), epsilon);

% Ensure dipolePerturbedRef is initialized before assigning values
dipolePerturbedRef = dip;
dipolePerturbedRef.complAmpl = normalizedPerturbedAmpReal + ...
                               1i * normalizedPerturbedAmpImag;

%% --- 1_e Validate Objective Function Before Optimization ---
ampGuess = dip.complAmpl; % dipolePerturbedRef.complAmpl;
dip.complAmpl = ampGuess; % Ensure ampGuess in dip
error = optimization.normObjectiveFunction_rad(dip, inputData);
disp(['Initial Test Error perturbed dip: ', num2str(error)]);
%% --- 2. Optimization Using PSO FF---
% --- Define Bounds ---
realLowerBounds = min(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
realUpperBounds = max(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
imagLowerBounds = min(normalizedPerturbedAmpImag) * ones(numDipoles, 1);
imagUpperBounds = max(normalizedPerturbedAmpImag) * ones(numDipoles, 1);

% Create final lower and upper bounds
lB = [realLowerBounds; imagLowerBounds];
uB = [realUpperBounds; imagUpperBounds];


% Create initial guess for PSO
initialGuess = [normalizedPerturbedAmpReal; normalizedPerturbedAmpImag]';

% Ensure correct size of initial swarm
swarmSize = numDipoles*2;
initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);

options_pso = optimoptions('particleswarm', ...
    'SwarmSize', swarmSize, ...                 % Number of particles in the swarm
    'MaxIterations', 200, ...                   % Maximum number of iterations
    'InertiaRange', [0.3, 1.5], ...             % Range for inertia weight (balance exploration/exploitation)
    'SelfAdjustmentWeight', 1.1, ...            % Weight for a particle's own best experience
    'SocialAdjustmentWeight', 1.05, ...         % Weight for following the global best solution
    'FunctionTolerance', 1e-10, ...             % Stop if function value improvement is below this threshold
    'MaxStallIterations', 40, ...               % Stop if no improvement in 40 consecutive iterations
    'InitialSwarmMatrix', initialSwarmMatrix,...% Set custom initial positions for the swarm
    'Display', 'iter' ...                       % Show progress at each iteration
    );  

optimFun = @(amp) ...
    optimization.normObjectiveFunction(amp(1:numDipoles).' + ...
    1i * amp((numDipoles+1):(numDipoles*2)).', dip, inputData);

% Run PSO for Amplitude Recovery
[optAmps_pso_vec, finalError_pso] = particleswarm(optimFun, 2 * numDipoles, lB, uB, options_pso);

optAmps_pso = nan(numDipoles,1);
optAmps_pso(:,1) = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec((numDipoles+1):(numDipoles*2));

disp(['Final Error (PSO): ', num2str(finalError_pso)]);

dipolePso = dipoleRef;
dipolePso.complAmpl = optAmps_pso;

fF_pso = fieldEvaluation.farField(rObserved, dipolePso, f0List);
totalPower_pso = sum(sum(fF_pso .* conj(fF_pso), 2) .* weights) / (2 * construct.Z0);

% PSO Plot
figure;
hold on;
plot(real(fF_pso/sqrt(totalPower_pso)), 'rx', 'MarkerSize', 6);
plot(imag(fF_pso)/sqrt(totalPower_pso), 'bx', 'MarkerSize', 6);
plot(real(fF_ref/sqrt(totalPower_ref)), 'ro', 'MarkerSize', 8);
plot(imag(fF_ref)/sqrt(totalPower_ref), 'bo', 'MarkerSize', 8);
title('PSO');
grid on;
hold off;





















% %% loading already generated fields
% nTh = size(farfield.theta, 2);
% nPh = size(farfield.phi, 1);
% 
% [x_theta, y_theta, z_theta] = utilities.transforms.sph2cartCoor(ones(nTh,1), farfield.theta.', ones(nTh,1) * farfield.phi(1));
% robs_the = [x_theta, y_theta, z_theta];
% 
% [x_phi, y_phi, z_phi] = utilities.transforms.sph2cartCoor(ones(nPh,1), ones(nPh,1) * farfield.theta(30), farfield.phi);
% robs_phi = [x_phi, y_phi, z_phi];
% 
% rObserved = [robs_the; robs_phi]; % Merge both planes
% 
% % theta = pi/2 to extract directivity in plane
% Directivity_phi_plane   = [farfield.D(30,:).', farfield.DTheta(30,:).', farfield.DPhi(30,:).'];
% Directivity_theta_plane = [farfield.D(:,1), farfield.DTheta(:,1), farfield.DPhi(:,1)];
% Directivity = [Directivity_theta_plane; Directivity_phi_plane];


% %% --- 1_c Compute Far-Field Parameters (Uniform Grid) ---
% % Define physical constants
% construct = utilities.constants.giveConstants();
% omega = 2 * pi * f0List;  % Angular frequency
% k = omega / construct.c0;  % Wavenumber
% rFar = 1e2 / k;            % Large observation distance
% 
% nTh = 60;  % Number of points for theta (0 to pi)
% nPh = 120;    % Number of points for phi (0 to 2pi)
% 
% % Create linearly spaced vectors
% theta = linspace(0, pi, nTh); 
% phi = linspace(0, 2*pi, nPh);
% 
% % Create a meshgrid
% [Phi, Theta] = meshgrid(phi,theta);
% 
% % Reshape into column vectors
% Theta = Theta(:); % Flatten into a single column vector
% Phi = Phi(:); % <=> Phi = reshape(Phi, [nTh*nPh, 1]);
% 
% % Convert spherical to Cartesian
% [x, y, z] = utilities.transforms.sph2cartCoor(ones(nTh*nPh, 1)*rFar, Theta, Phi);
% rObserved = [x, y, z]; % Combine Cartesian coordinates into a matrix
% 
% % Evaluate the far field
% [eF] = fieldEvaluation.eleFieldM2(rObserved, dip, f0List);
% [mF] = fieldEvaluation.magFieldM2(rObserved, dip, f0List);
% [fF] = fieldEvaluation.farField(rObserved, dip, f0List);
% [propFact] = fieldEvaluation.computePropagationFactor(rObserved,k);
% eFfF = fF.*propFact;
% 
% power = fieldEvaluation.powerPoynting(eF,mF);
% powerFar = fieldEvaluation.powerDensityFar(eFfF,rObserved,dip);
% [fF2(:, 1), fF2(:, 2), fF2(:, 3)] = utilities.transforms.sph2cartVec(Theta,Phi,ones(nTh*nPh,1),farfield.FTheta(:),farfield.FPhi(:));
% 
% % Transform Cartesian field components to spherical
% [~, Fth, Fph] = utilities.transforms.cart2sphVec(...
%                 x, y, z, fF(:, 1), fF(:, 2), fF(:, 3) ...
%                                                  );
% 
% % Reshape Theta, Phi, Fth, and Fph for contour plotting
% Theta = reshape(Theta, [nTh, nPh]);
% Phi = reshape(Phi, [nTh, nPh]);
% Fth = reshape(Fth, [nTh, nPh]);
% Fph = reshape(Fph, [nTh, nPh]);
% 
% % Visualization of Fth and Fph using contour plots with color bars
% 
% % Plot for the real part of Fth
% figure
% contourf(Theta/pi, Phi/(pi), real(Fth));
% xlabel('\theta/ \pi')
% ylabel('\phi/ \pi')
% title('Re[Fth]');
% colorbar; % Add color bar for reference
% 
% 
% % Plot for the imaginary part of Fth
% figure
% contourf(Theta/(pi), Phi/pi, imag(Fth));
% xlabel('\theta/ \pi')
% ylabel('\phi/ \pi')
% title('Im[Fth]');
% colorbar;
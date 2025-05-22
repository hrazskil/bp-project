%% --- 1_a Initialization ---
clc; clear; close all;

% Load dipole variables from file 
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\halfwaveDipole.mat');
frequency = f0List; % Frequency of 1 GHz

dipoleRef = dip;

% Extract numDipoles
numDipoles = numel(dip.complAmpl);

% Compute the maximum magnitude across all dipole amplitudes
% maxAmp = max(abs(dip.complAmpl));
% Normalize the complex amplitudes
% dipoleRef.complAmpl = dip.complAmpl / maxAmp;

%% Perturb amplitudes
% Generate complex perturbation factors with small variations
% perturbationFactor = 1 + 0.1 *(randn(numDipoles, 1));
% Apply perturbations to the normalized complex amplitudes
% perturbedAmp = dipoleRef.complAmpl .* perturbationFactor;
% Normalize the perturbed amplitudes
% maxPerturbedAmp = max(abs(perturbedAmp));
% dipolePerturbed.complAmpl = perturbedAmp / maxPerturbedAmp;

perturbationFactor = 0.75*abs(dipoleRef.complAmpl).*(randn(numDipoles, 1)+1i*randn(numDipoles, 1));

dipolePerturbed = dipoleRef;

dipolePerturbed.complAmpl = dipoleRef.complAmpl + perturbationFactor;
%% --- 1_c Compute Far-Field Parameters ---
% Define physical constants
construct = utilities.constants.giveConstants();
omega = 2 * pi * frequency;  % Angular frequency
k = omega / construct.c0;    % Wavenumber
rFar = 1e6 / k;              % Large observation distance

% degree: { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,
%     350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,
%     3470, 3890, 4334, 4802, 5294, 5810 };
Nleb = 350;                 % Number of Lebedev quadrature points

% Get Lebedev quadrature points and weights
[points, weights, ~] = utilities.getLebedevSphere(Nleb);
rObserved = points * rFar;   % Scale points to observation distance

% Compute far-field patterns
fF_ref = fieldEvaluation.farFieldM2(rObserved, dipoleRef, frequency);
fF_Perturbed = fieldEvaluation.farFieldM2(rObserved, dipolePerturbed, frequency);

% Compute total radiated power for normalization
totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);
totalPower_Perturbed = sum(sum(fF_Perturbed .* conj(fF_Perturbed), 2) .* weights) / (2 * construct.Z0);

% normalize
dipoleRef.complAmpl = dipoleRef.complAmpl / sqrt(totalPower_ref);
dipolePerturbed.complAmpl = dipolePerturbed.complAmpl / sqrt(totalPower_Perturbed);

%% === Far-Field Comparison: Perturbed vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipolePerturbed, frequency, 180, 360);

%% === Far-Field Comparison: Reference vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleRef, frequency, 180, 360);

%% === Far-Field Comparison: Perturbed vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipolePerturbed, frequency, 180, 360);

%% === Far-Field Intensity Comparison: Perturbed vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipolePerturbed, frequency, 180, 360);

%% --- Validate Objective Function Before Optimization ---
initialError = optimization.normObjectiveFunction( ...
    dipolePerturbed.complAmpl, dip, frequency, points,  weights, ...
    fF_ref, totalPower_ref);
disp(['Initial Test Error: ', num2str(initialError)]);

% --- Define Bounds ---
realPert = real(dipolePerturbed.complAmpl);
imagPert = imag(dipolePerturbed.complAmpl);

ampMinReal = min(realPert);  ampMaxReal = max(realPert);
ampMinImag = min(imagPert);  ampMaxImag = max(imagPert);

lB = [ampMinReal * ones(numDipoles, 1); ampMinImag * ones(numDipoles, 1)];
uB = [ampMaxReal * ones(numDipoles, 1); ampMaxImag * ones(numDipoles, 1)];


%% --- 2. Optimization Using PSO ---
% Create initial guess for PSO
initialGuess = [real(dipolePerturbed.complAmpl); imag(dipolePerturbed.complAmpl)]';

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
    1i * amp((numDipoles+1):(numDipoles*2)).', dip, frequency, points, ...
    weights, fF_ref, totalPower_ref);

% Run PSO for Amplitude Recovery
[optAmps_pso_vec, finalError_pso] = particleswarm(optimFun, 2 * numDipoles, lB, uB, options_pso);

optAmps_pso = nan(numDipoles,1);
optAmps_pso(:,1) = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec((numDipoles+1):(numDipoles*2));

disp(['Final Error (PSO): ', num2str(finalError_pso)]);

dipolePso = dipoleRef;
dipolePso.complAmpl = optAmps_pso;

%% --- 3. Optimization Using fmincon ---

initialGuess = [optAmps_pso_vec(1:numDipoles);...                          % serialization of optimizations
              optAmps_pso_vec(numDipoles+1:end)]';

options_fmincon = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...             % Use Sequential Quadratic Programming (SQP) algorithm
    'MaxIterations', 100, ...           % Maximum number of iterations
    'MaxFunctionEvaluations', 5000, ... % Maximum number of function evaluations
    'OptimalityTolerance', 1e-10, ...   % Stop if optimality conditions are met within this tolerance
    'StepTolerance', 1e-10, ...         % Stop if step size is below this threshold
    'Display', 'iter');                 % Display iteration details 

optimFun_fmincon = @(amp) ...
    optimization.normObjectiveFunction(amp(1:numDipoles).' + ...
    1i * amp(numDipoles+1:end).', dip, frequency, points, weights, ...
    fF_ref, totalPower_ref);

[optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun_fmincon,...
    initialGuess, [], [], [], [], lB, uB, [], options_fmincon);

optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                  1i * optAmps_fmincon_vec(numDipoles+1:end).';
disp(['Final Error (fmincon): ', num2str(finalError_fmincon)]);


dipoleFmincon = dipoleRef;
dipoleFmincon.complAmpl = optAmps_fmincon;

fF_Fmincon = fieldEvaluation.farFieldM2(rObserved, dipoleFmincon, frequency);
totalPower_Fmincon = sum(sum(fF_Fmincon .* conj(fF_Fmincon), 2) .* weights) / (2 * construct.Z0);






% Fmincon Plot
%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);


%% Step 1: Define Angular Grid
theta = linspace(0, pi, Ntheta).';              % Zenith angle (0 to pi)
phi   = linspace(0, 2*pi, Nphi);                % Azimuth angle (0 to 2pi)
[PHI, THETA] = meshgrid(phi, theta);            % Create 2D angular grid

%% Step 2: Convert to Cartesian Unit Vectors
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta * Nphi, 1), THETA(:), PHI(:));
r_hat = [x, y, z];                    % Flatten to N x 3 matrix

%% Step 3: Evaluate Far-Field Electric Field Vectors
fF_ref = fieldEvaluation.farFieldM2(r_hat, dipoleRef, freq);
fF_pso = fieldEvaluation.farFieldM2(r_hat, dipolePso, freq);

%% Step 4: Convert Cartesian to Spherical Field Components
[~, Fth_ref, Fph_ref] = utilities.transforms.cart2sphVec(x(:), y(:), z(:), fF_ref(:,1), fF_ref(:,2), fF_ref(:,3));
[~, Fth_pso, Fph_pso] = utilities.transforms.cart2sphVec(x(:), y(:), z(:), fF_pso(:,1), fF_pso(:,2), fF_pso(:,3));

%% Step 5: Reshape Components to Angular Grid
Fth_ref = reshape(Fth_ref, Ntheta, Nphi);
Fph_ref = reshape(Fph_ref, Ntheta, Nphi);
Fth_pso = reshape(Fth_pso, Ntheta, Nphi);
Fph_pso = reshape(Fph_pso, Ntheta, Nphi);

%% Step 6: Compute Total Intensity
intensity_ref = abs(Fth_ref).^2 + abs(Fph_ref).^2;
intensity_pso = abs(Fth_pso).^2 + abs(Fph_pso).^2;

%% Step 7: Normalize (Linear Scale)
intensity_ref_norm = intensity_ref / max(intensity_ref(:));
intensity_pso_norm = intensity_pso / max(intensity_pso(:));
intensity_diff = abs(intensity_pso_norm - intensity_ref_norm);
%% --- 4. Near-Field Evaluation and Visualization ---

% Define a 2D grid in the XZ-plane centered around the dipole structure
gridRange = 0.2;            % Spatial extent in meters (adjust as needed)
gridRes = 200;              % Number of grid points per axis

xVals = linspace(-gridRange, gridRange, gridRes);
yVals = 0;                  % Single plane: y = 0 (XZ plane)
zVals = linspace(-gridRange, gridRange, gridRes);

gridSpec.x = xVals;
gridSpec.y = yVals;
gridSpec.z = zVals;

% Evaluate near-field using optimized dipole model
[Efield, Hfield, Sfield, grid] = utilities.visualizations.evaluateNearFieldGrid(dipoleFmincon, frequency, gridSpec);

% Visualize Poynting vector magnitude (log scale recommended for clarity)
figure;
imagesc(grid.x, grid.z, log10(squeeze(Sfield(:,1,:))));
axis equal tight;
xlabel('x [m]');
ylabel('z [m]');
title('Near-Field Poynting Vector Magnitude (log_{10})');
colorbar;
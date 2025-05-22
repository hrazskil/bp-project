%% --- 1_b Initialization ---
clc; clear; close all;

%% Load precomputed far-field data and dipole configuration
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\dataVertical.tsv')
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\dataHorizontal.tsv')
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\halfwaveDipole.mat')

% Determine number of observation points from data files
Nbeta = size(dataHorizontal(:,2),1);  
Nphi = size(dataVertical(:,2),1);    

% Reconstruct observation point vectors from saved angular values
x1 = cos(dataVertical(:,1));
y1 = 0*x1;
z1 = sin(dataVertical(:,1));

x2 = cos(dataHorizontal(:,1));
y2 = sin(dataHorizontal(:,1));
z2 = 0*x2;

vertical_points = [x1(:,1), y1(:,1), z1(:,1)];
horizontal_points = [x2(:,1), y2(:,1), z2(:,1)];

% Assemble structured input data for optimization
inputData.horizontal.rad        = dataHorizontal(:,2);
inputData.horizontal.points     = horizontal_points;
inputData.horizontal.weights    = ones(Nbeta, 1);  % Uniform weights

inputData.horizontal.totalPower = ...
    sum((inputData.horizontal.rad .* inputData.horizontal.weights));

inputData.vertical.rad          = dataVertical(:,2);
inputData.vertical.points       = vertical_points;
inputData.vertical.weights      = ones(Nbeta, 1);  % Uniform weights

inputData.vertical.totalPower = ...
    sum((inputData.vertical.rad .* inputData.vertical.weights));

inputData.freq = f0List;  % Frequency list used in computation

%% Reference dipole amplitudes
dipoleRef = dip;

% Extract numDipoles
numDipoles = numel(dipoleRef.complAmpl);

%% Compute initial error with true amplitudes
ampGuess = dip.complAmpl; 
dip.complAmpl = ampGuess; 
error = optimization.normObjectiveFunction_rad(dipoleRef, inputData);
disp(['Initial Test After Normalization: ', num2str(error)]);

%% Perturb amplitudes
perturbationFactor = 0.5*abs(dipoleRef.complAmpl).*(randn(numDipoles, 1)+1i*randn(numDipoles, 1));

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

%% --- 2. Optimization Using PSO FF---

% Use perturbed amplitudes for error evaluation
dip.complAmpl = dipolePerturbed.complAmpl;
error = optimization.normObjectiveFunction_rad(dip, inputData);
disp(['Initial Test Error perturbed dip: ', num2str(error)]);

% --- Define Bounds ---
realPert = real(dipolePerturbed.complAmpl);
imagPert = imag(dipolePerturbed.complAmpl);

ampMinReal = min(realPert);  ampMaxReal = max(realPert);
ampMinImag = min(imagPert);  ampMaxImag = max(imagPert);

lB = [ampMinReal * ones(numDipoles, 1); ampMinImag * ones(numDipoles, 1)];
uB = [ampMaxReal * ones(numDipoles, 1); ampMaxImag * ones(numDipoles, 1)];

% --- Initialize PSO Parameters ---
% Combine real and imaginary parts into a single initial guess vector
% Create initial guess for PSO
initialGuess = [real(dipolePerturbed.complAmpl); imag(dipolePerturbed.complAmpl)]';

% Create swarm matrix by replicating the initial guess (each row is a particle)
swarmSize = numDipoles * 2;  % Swarm size is double the number of dipoles (real + imag)
initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);

% Define PSO optimization settings
options_pso = optimoptions('particleswarm', ...
    'SwarmSize', swarmSize, ...                 % Number of particles
    'MaxIterations', 200, ...                   % Iteration limit
    'InertiaRange', [0.3, 1.5], ...             % Inertia control for convergence behavior
    'SelfAdjustmentWeight', 1.1, ...            % Particle's self-exploration factor
    'SocialAdjustmentWeight', 1.05, ...         % Attraction to global best solution
    'FunctionTolerance', 1e-10, ...             % Stop when improvement is below threshold
    'MaxStallIterations', 40, ...               % Stop when no progress in 40 iterations
    'InitialSwarmMatrix', initialSwarmMatrix,...% Custom initial particle positions
    'Display', 'iter' ...                       % Show iteration info in console
    );  

% --- Objective Function Definition ---
% This function evaluates the fitness of a candidate amplitude vector
% Split real and imag parts, construct complex amplitude vector, and evaluate error
optimFun = @(amp) optimization.optimFunX(amp, dip, inputData, numDipoles);

% --- Run Particle Swarm Optimization ---
% Optimize dipole amplitudes to match far-field radiation
[optAmps_pso_vec, finalError_pso] = particleswarm(optimFun, 2 * numDipoles, lB, uB, options_pso);

% Combine optimized real and imaginary parts into complex amplitudes
optAmps_pso = nan(numDipoles,1);
optAmps_pso(:,1) = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec((numDipoles+1):(numDipoles*2));

% Show final optimization error
disp(['Final Error (PSO): ', num2str(finalError_pso)]);
dipolePso = dipoleRef;
dipolePso.complAmpl=optAmps_pso;


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

optimFun = @(amp) optimization.optimFunX(amp, dip, inputData, numDipoles);

[optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun,...
    initialGuess, [], [], [], [], lB, uB, [], options_fmincon);

optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                  1i * optAmps_fmincon_vec(numDipoles+1:end).';
disp(['Final Error (fmincon): ', num2str(finalError_fmincon)]);


dipoleFmincon1 = dipoleRef;
dipoleFmincon1.complAmpl = optAmps_fmincon;

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon1, inputData.freq, 180, 360);


%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon1, inputData.freq, 180, 360);


% old nearfield
% %% Parameters for grid
% nTh = 181; nPh = 361;
% theta = linspace(0, pi, nTh); 
% phi = linspace(0, 2*pi, nPh);
% [Phi, Theta] = meshgrid(phi, theta);
% x = sin(Theta) .* cos(Phi);
% y = sin(Theta) .* sin(Phi);
% z = cos(Theta);
% rObs = [x(:), y(:), z(:)];
% 
% % Compute fields
% E_ref = fieldEvaluation.farFieldM2(rObs, dipoleRef, f0List);
% E_pso = fieldEvaluation.farFieldM2(rObs, dipolePso, f0List);
% E_fmc = fieldEvaluation.farFieldM2(rObs, dipoleFmincon, f0List);
% 
% % Convert to power density
% I_ref = sum(abs(E_ref).^2, 2);
% I_pso = sum(abs(E_pso).^2, 2);
% I_fmc = sum(abs(E_fmc).^2, 2);
% 
% % Normalize
% I_ref = I_ref / max(I_ref);
% I_pso = I_pso / max(I_pso);
% I_fmc = I_fmc / max(I_fmc);
% 
% % Absolute difference
% diff_pso = abs(I_pso - I_ref);
% diff_fmc = abs(I_fmc - I_ref);
% 
% % Reshape to 2D for plotting
% diff_pso_grid = reshape(diff_pso, [nTh, nPh]);
% diff_fmc_grid = reshape(diff_fmc, [nTh, nPh]);
% 
% % Plot PSO error
% figure;
% contourf(theta/pi, phi/pi, diff_pso_grid, 50, 'LineColor','none');
% colorbar;
% xlabel('\theta / \pi'); ylabel('\phi / \pi');
% title('PSO Optimization – Intensity Difference (|I_{opt} - I_{ref}|)');
% 
% % Plot fmincon error
% figure;
% contourf(theta/pi, phi/pi, diff_fmc_grid, 50, 'LineColor','none');
% colorbar;
% xlabel('\theta / \pi'); ylabel('\phi / \pi');
% title('fmincon Optimization – Intensity Difference (|I_{opt} - I_{ref}|)');
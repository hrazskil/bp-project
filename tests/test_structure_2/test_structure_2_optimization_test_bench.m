%% --- Initialization for test_structure_2 ---
clc; clear; close all;

%% Load saved far-field data and dipole configuration
load('tests/test_structure_2/dataVertical.tsv');   % vertical far-field data
load('tests/test_structure_2/dataHorizontal.tsv'); % horizontal far-field data
load('tests/test_structure_2/DipoleArray.mat');    % dipole geometry + f0List

Nbeta = size(dataHorizontal(:,2),1);  % number of azimuthal observations
Nphi = size(dataVertical(:,2),1);     % number of elevation observations

% Reconstruct Cartesian observation points from angles
x1 = cos(dataVertical(:,1));
y1 = 0*x1;
z1 = sin(dataVertical(:,1));

x2 = cos(dataHorizontal(:,1));
y2 = sin(dataHorizontal(:,1));
z2 = 0*x2;

vertical_points = [x1, y1, z1];
horizontal_points = [x2, y2, z2];

% Create input data structure for optimization
inputData.horizontal.rad = dataHorizontal(:,2);
inputData.horizontal.points = horizontal_points;
inputData.horizontal.weights = ones(Nbeta, 1);
inputData.horizontal.totalPower = sum(inputData.horizontal.rad);

inputData.vertical.rad = dataVertical(:,2);
inputData.vertical.points = vertical_points;
inputData.vertical.weights = ones(Nphi, 1);
inputData.vertical.totalPower = sum(inputData.vertical.rad);

inputData.freq = f0List;

%% Compute initial error with true amplitudes
ampGuess = dip.complAmpl; 
dip.complAmpl = ampGuess; 
tic
error = optimization.normObjectiveFunction_rad(dip, inputData);
toc
disp(['Initial Test Error: ', num2str(error)]);

%% Normalize true dipole amplitudes
numDipoles = numel(dip.complAmpl);
dipoleAmpReal = real(dip.complAmpl);
dipoleAmpImag = imag(dip.complAmpl);

maxAmpReal = max(max(abs(dipoleAmpReal)),eps);
maxAmpImag = max(max(abs(dipoleAmpImag)),eps);

normalizedAmpReal = dipoleAmpReal / maxAmpReal;
normalizedAmpImag = dipoleAmpImag / maxAmpReal;

dipoleRef = dip;
dipoleRef.complAmpl = normalizedAmpReal + 1i * normalizedAmpImag;

tic
error = optimization.normObjectiveFunction_rad(dipoleRef, inputData);
toc
disp(['Initial Test After Normalization: ', num2str(error)]);

%% Create perturbed (test) amplitudes with small noise
realPerturbationFactor = 1 + 0.1 * randn(numDipoles, 1);
imagPerturbationFactor = 1 + 0.1 * randn(numDipoles, 1);

perturbedAmp = real(dip.complAmpl) .* realPerturbationFactor + ...
               1i * imag(dip.complAmpl) .* imagPerturbationFactor;

% Normalize perturbed amplitudes
perturbedAmpReal = real(perturbedAmp);
perturbedAmpImag = imag(perturbedAmp);
normalizedPerturbedAmpReal = perturbedAmpReal ./ max(max(abs(perturbedAmpReal)),eps);
normalizedPerturbedAmpImag = perturbedAmpImag ./ max(max(abs(perturbedAmpReal)),eps);

% Create new dipole object with perturbed amplitudes
dipolePerturbedRef = dip;
dipolePerturbedRef.complAmpl = normalizedPerturbedAmpReal + ...
                               1i * normalizedPerturbedAmpImag;

% Use perturbed amplitudes for error evaluation
dip.complAmpl = dipolePerturbedRef.complAmpl;
tic
error = optimization.normObjectiveFunction_rad(dip, inputData);
toc
disp(['Initial Test Error perturbed dip: ', num2str(error)]);

%% === Far-Field Comparison: Perturbed vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipolePerturbedRef, inputData.freq, 180, 360);

%% === Far-Field Intensity Comparison: Perturbed vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipolePerturbedRef, inputData.freq, 180, 360, [2 98], [2 98], [1 99], 0.0005);
%% --- 2. Optimization Using PSO FF---

% --- Define Bounds ---
% Define element-wise lower and upper bounds for the real and imaginary parts 
% of the dipole amplitudes based on normalized perturbed values
realLowerBounds = min(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
realUpperBounds = max(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
imagLowerBounds = min(normalizedPerturbedAmpImag) * ones(numDipoles, 1);
imagUpperBounds = max(normalizedPerturbedAmpImag) * ones(numDipoles, 1);

% Combine bounds into single vectors for PSO input
lB = [realLowerBounds; imagLowerBounds];
uB = [realUpperBounds; imagUpperBounds];

% --- Initialize PSO Parameters ---
% Combine real and imaginary parts into a single initial guess vector
initialGuess = [normalizedPerturbedAmpReal; normalizedPerturbedAmpImag]';

% Create swarm matrix by replicating the initial guess (each row is a particle)
swarmSize = numDipoles/4;  % Swarm size is double the number of dipoles (real + imag)
initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);

% Define PSO optimization settings
options_pso = optimoptions('particleswarm', ...
    'SwarmSize', swarmSize, ...                 % Number of particles
    'MaxIterations', 400, ...                   % Iteration limit
    'InertiaRange', [0.3, 1.5], ...             % Inertia control for convergence behavior
    'SelfAdjustmentWeight', 1.1, ...            % Particle's self-exploration factor
    'SocialAdjustmentWeight', 1.05, ...         % Attraction to global best solution
    'FunctionTolerance', 1e-6, ...             % Stop when improvement is below threshold
    'MaxStallIterations', 20, ...               % Stop when no progress in 40 iterations
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

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipolePso, inputData.freq, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipolePso, inputData.freq, 180, 360, [2 98], [2 98], [1 99], 0.0005);

%% --- 3. Optimization Using fmincon ---
initialGuess = [optAmps_pso_vec(1:numDipoles);...                          % serialization of optimizations
              optAmps_pso_vec(numDipoles+1:end)]';

options_fmincon = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...             % Use Sequential Quadratic Programming (SQP) algorithm
    'MaxIterations', 100, ...           % Maximum number of iterations
    'MaxFunctionEvaluations', 50000, ... % Maximum number of function evaluations
    'OptimalityTolerance', 1e-12, ...   % Stop if optimality conditions are met within this tolerance
    'StepTolerance', 1e-12, ...         % Stop if step size is below this threshold
    'Display', 'iter');                 % Display iteration details 

optimFun = @(amp) optimization.optimFunX(amp, dip, inputData, numDipoles);

[optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun,...
    initialGuess, [], [], [], [], lB, uB, [], options_fmincon);

optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                  1i * optAmps_fmincon_vec(numDipoles+1:end).';
disp(['Final Error (fmincon): ', num2str(finalError_fmincon)]);


dipoleFmincon = dipoleRef;
dipoleFmincon.complAmpl = optAmps_fmincon;

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360);


%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360, [2 98], [2 98], [1 99], 0.00005);



%% PSO optimization
options_pso = optimoptions('particleswarm', ...
    'SwarmSize', 100, 'MaxIterations', 200, 'Display', 'iter');

[x_pso, err_pso] = particleswarm(optimFun, 2*numDipoles, lb, ub, options_pso);
disp(['PSO error: ', num2str(err_pso)]);

%% fmincon refinement
options_fmc = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'MaxIterations', 300, 'MaxFunctionEvaluations', 5000, 'Display', 'iter');

[x_fmc, err_fmc] = fmincon(optimFun, x_pso, [], [], [], [], lb, ub, [], options_fmc);
disp(['Final error after fmincon: ', num2str(err_fmc)]);

% Reconstruct final complex amplitudes
dipOptimized = dip;
dipOptimized.complAmpl = x_fmc(1:numDipoles) * maxReal + 1i * x_fmc(numDipoles+1:end) * maxImag;

%% Final evaluation
finalError = optimization.normObjectiveFunction_rad(dipOptimized, inputData);
disp(['Final normalized error: ', num2str(finalError)]);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360);


%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360, [2 98], [2 98], [1 99], 0.00005);

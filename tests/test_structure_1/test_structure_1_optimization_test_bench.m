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

%% Compute initial error with true amplitudes
ampGuess = dip.complAmpl; 
dip.complAmpl = ampGuess; 
error = optimization.normObjectiveFunction_rad(dip, inputData);
disp(['Initial Test Error: ', num2str(error)]);

%% Normalize true dipole amplitudes
numDipoles = numel(dip.complAmpl);
dipoleAmpReal = real(dip.complAmpl);
dipoleAmpImag = imag(dip.complAmpl);

maxAmpReal = max(max(abs(dipoleAmpReal)),eps);
maxAmpImag = max(max(abs(dipoleAmpImag)),eps);

normalizedAmpReal = dipoleAmpReal / maxAmpReal;
normalizedAmpImag = dipoleAmpImag / maxAmpImag;

dipoleRef = dip;
dipoleRef.complAmpl = normalizedAmpReal + 1i * normalizedAmpImag;

%% Compute initial error with true amplitudes
ampGuess = dip.complAmpl; 
dip.complAmpl = ampGuess; 
error = optimization.normObjectiveFunction_rad(dipoleRef, inputData);
disp(['Initial Test After Normalization: ', num2str(error)]);

%% Create perturbed (test) amplitudes with small noise
realPerturbationFactor = 1 + 0.01 * randn(numDipoles, 1);
imagPerturbationFactor = 1 + 0.01 * randn(numDipoles, 1);

perturbedAmp = real(dip.complAmpl) .* realPerturbationFactor + ...
               1i * imag(dip.complAmpl) .* imagPerturbationFactor;

% Normalize perturbed amplitudes
perturbedAmpReal = real(perturbedAmp);
perturbedAmpImag = imag(perturbedAmp);

normalizedPerturbedAmpReal = perturbedAmpReal ./ max(max(abs(perturbedAmpReal)),eps);
normalizedPerturbedAmpImag = perturbedAmpImag ./ max(max(abs(perturbedAmpImag)),eps);

% Create new dipole object with perturbed amplitudes
dipolePerturbedRef = dip;
dipolePerturbedRef.complAmpl = normalizedPerturbedAmpReal + ...
                               1i * normalizedPerturbedAmpImag;

% Use perturbed amplitudes for error evaluation
dip.complAmpl = dipolePerturbedRef.complAmpl;
error = optimization.normObjectiveFunction_rad(dip, inputData);
disp(['Initial Test Error perturbed dip: ', num2str(error)]);

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

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipolePso, inputData.freq, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipolePso, inputData.freq, 180, 360, [2 98], [2 98], [1 99], 0.0005);

%% --- 3. Optimization Using fmincon ---
initialGuess = [optAmps_pso_vec(1:numDipoles);...                          % serialization of optimizations
              optAmps_pso_vec(numDipoles+1:end)]';

options_fmincon = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...             % Use Sequential Quadratic Programming (SQP) algorithm
    'MaxIterations', 100, ...           % Maximum number of iterations
    'MaxFunctionEvaluations', 5000, ... % Maximum number of function evaluations
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
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360, [2 98], [2 98], [1 99], 0.00005);
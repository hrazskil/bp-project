%% --- 1_a Initialization ---
clc; clear; close all;

% Load dipole variables from file 
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\halfwaveDipole.mat');
frequency = f0List; % Frequency of 1 GHz

% Extract real and imaginary parts of dipole amplitudes
numDipoles = numel(dip.complAmpl);
dipoleAmpReal = real(dip.complAmpl);
dipoleAmpImag = imag(dip.complAmpl);

% Normalize the real and imaginary parts separately and keep norms
maxAmpReal = max(abs(dipoleAmpReal));
maxAmpImag = max(abs(dipoleAmpImag));
normalizedAmpReal = dipoleAmpReal / maxAmpReal;
normalizedAmpImag = dipoleAmpImag / maxAmpImag;

% Save reference dipole amplitudes after normalization
dipoleRef = dip;
dipoleRef.complAmpl = normalizedAmpReal + 1i * normalizedAmpImag;


%% --- 1_b Perturbation of Initial Amplitudes ---
% Generate independent perturbation factors
realPerturbationFactor = 1 + 0.1 * randn(numDipoles, 1);
imagPerturbationFactor = 1 + 0.1 * randn(numDipoles, 1);

% Apply perturbations separately to real and imaginary parts
perturbedAmp = real(dip.complAmpl) .* realPerturbationFactor + ...
               1i * imag(dip.complAmpl) .* imagPerturbationFactor;

% Normalize perturbed dipole amplitudes
perturbedAmpReal = real(perturbedAmp);
perturbedAmpImag = imag(perturbedAmp);
normalizedPerturbedAmpReal = perturbedAmpReal / max(abs(perturbedAmpReal));
normalizedPerturbedAmpImag = perturbedAmpImag / max(abs(perturbedAmpImag));

dipolePerturbedRef.complAmpl = normalizedPerturbedAmpReal + ...
                               1i * normalizedPerturbedAmpImag;

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

% Compute reference far-field pattern
fF_ref = fieldEvaluation.farFieldM2(rObserved, dipoleRef, frequency);

% Compute total radiated power for normalization
totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);


%% --- 1_d Validate Objective Function Before Optimization ---
initialError = optimization.normObjectiveFunction( ...
    dipolePerturbedRef.complAmpl, dip, frequency, points, weights, ...
    fF_ref, totalPower_ref);
disp(['Initial Test Error: ', num2str(initialError)]);

% --- Define Bounds ---
realLowerBounds = min(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
realUpperBounds = max(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
imagLowerBounds = min(normalizedPerturbedAmpImag) * ones(numDipoles, 1);
imagUpperBounds = max(normalizedPerturbedAmpImag) * ones(numDipoles, 1);

% Create final lower and upper bounds
lB = [realLowerBounds; imagLowerBounds];
uB = [realUpperBounds; imagUpperBounds];


%% --- 2. Optimization Using PSO ---
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
    1i * amp((numDipoles+1):(numDipoles*2)).', dip, frequency, points, ...
    weights, fF_ref, totalPower_ref);

% Run PSO for Amplitude Recovery
[optAmps_pso_vec, finalError_pso] = particleswarm(optimFun, 2 * numDipoles, lB, uB, options_pso);

optAmps_pso = nan(numDipoles,1);
optAmps_pso(:,1) = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec((numDipoles+1):(numDipoles*2));

disp(['Final Error (PSO): ', num2str(finalError_pso)]);

dipolePso = dipoleRef;
dipolePso.complAmpl = optAmps_pso;

fF_pso = fieldEvaluation.farFieldM2(rObserved, dipolePso, frequency);
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

%% --- 3. Optimization Using fmincon ---

% initialGuess = [optAmps_pso_vec(1:numDipoles);...                          % serialization of optimizations
%               optAmps_pso_vec(numDipoles+1:end)]';

options_fmincon = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...             % Use Sequential Quadratic Programming (SQP) algorithm
    'MaxIterations', 100, ...           % Maximum number of iterations
    'MaxFunctionEvaluations', 5000, ... % Maximum number of function evaluations
    'OptimalityTolerance', 1e-12, ...   % Stop if optimality conditions are met within this tolerance
    'StepTolerance', 1e-12, ...         % Stop if step size is below this threshold
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
figure;
hold on;
plot(real(fF_Fmincon/sqrt(totalPower_Fmincon)), 'rx', 'MarkerSize', 6);
plot(imag(fF_Fmincon)/sqrt(totalPower_Fmincon), 'bx', 'MarkerSize', 6);
plot(real(fF_ref/sqrt(totalPower_ref)), 'ro', 'MarkerSize', 8);
plot(imag(fF_ref)/sqrt(totalPower_ref), 'bo', 'MarkerSize', 8);
title('Fmincon');
grid on;
hold off;

%% --- 4. Optimization Using GA ---

initialGuess = [normalizedPerturbedAmpReal; normalizedPerturbedAmpImag]';

% Ensure correct size of initial population
populationSize = swarmSize;
initialPopulationMatrix = repmat(initialGuess, populationSize, 1);

options_ga = optimoptions('ga', ...
    'PopulationSize', populationSize, ...      % Number of individuals in the population
    'MaxGenerations', 200, ...                 % Maximum number of generations
    'CrossoverFraction', 0.8, ...              % Crossover fraction
    'EliteCount', 2, ...                       % Number of elite individuals
    'FunctionTolerance', 1e-8, ...             % Stop if function value improvement is below this threshold
    'MutationFcn', {@mutationadaptfeasible, 0.1}, ... % Mutation function
    'Display', 'iter', ...                     % Show progress at each iteration
    'InitialPopulationMatrix', initialPopulationMatrix ... % Set custom initial population
    );

optimFun_ga = @(amp) ...
    optimization.normObjectiveFunction(amp(1:numDipoles).' + ...
    1i * amp((numDipoles+1):(numDipoles*2)).', dip, frequency, points, ...
    weights, fF_ref, totalPower_ref);

% Run GA for Amplitude Recovery
[optAmps_ga_vec, finalError_ga] = ga(optimFun_ga, 2 * numDipoles, [], [], [], [], lB, uB, [], options_ga);

optAmps_ga = nan(numDipoles, 1);
optAmps_ga(:, 1) = optAmps_ga_vec(1:numDipoles) + 1i * optAmps_ga_vec((numDipoles+1):(numDipoles*2));

disp(['Final Error (GA): ', num2str(finalError_ga)]);

dipoleGa = dipoleRef;
dipoleGa.complAmpl = optAmps_ga;

fF_ga = fieldEvaluation.farFieldM2(rObserved, dipoleGa, frequency);
totalPower_ga = sum(sum(fF_ga .* conj(fF_ga), 2) .* weights) / (2 * construct.Z0);

% GA Plot
figure;
hold on;
plot(real(fF_ga/sqrt(totalPower_ga)), 'rx', 'MarkerSize', 6);
plot(imag(fF_ga)/sqrt(totalPower_ga), 'bx', 'MarkerSize', 6);
plot(real(fF_ref/sqrt(totalPower_ref)), 'ro', 'MarkerSize', 8);
plot(imag(fF_ref)/sqrt(totalPower_ref), 'bo', 'MarkerSize', 8);
title('GA');
grid on;
hold off;
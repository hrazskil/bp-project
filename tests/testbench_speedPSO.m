%% --- PSO Optimization for Different Numbers of Dipoles ---
clc; clear; close all;

% Load dipole variables from file
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat')

% Define available number of dipoles for testing
dipoleCounts = round(linspace(1, numel(dip.complAmpl)/4, 4));  % Varying dipole counts (from 10 to the total number of dipoles)

% Define fixed Lebedev quadrature degree (e.g., degree 302 for now)
degree = 302;

% Define physical constants
construct = utilities.constants.giveConstants();
omega = 2 * pi * f0List;  % Angular frequency
k = omega / construct.c0;    % Wavenumber
rFar = 1e6 / k;              % Large observation distance

% Generate Lebedev quadrature points and weights for fixed degree
[points, weights, ~] = utilities.getLebedevSphere(degree);
rObserved = points * rFar;  % Scale points to observation distance

% Compute reference far-field pattern for fixed Lebedev 302
fF_ref = fieldEvaluation.farFieldM2(rObserved, dip, f0List);
totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);

% Time measurement for all dipole subsets
timeResults = zeros(length(dipoleCounts), 1);
finalErrorResults = zeros(length(dipoleCounts), 1);

for i = 1:length(dipoleCounts)
    % Subset dipoles based on current count
    numDipoles = dipoleCounts(i);
    dipSubset.complAmpl = dip.complAmpl(1:numDipoles);
    dipSubset.pos = dip.pos(1:numDipoles, :);
    dipSubset.dir = dip.dir(1:numDipoles, :);
    
    % Perturbation of Initial Amplitudes
    realPerturbationFactor = 1 + 0.1 * randn(numDipoles, 1);
    imagPerturbationFactor = 1 + 0.1 * randn(numDipoles, 1);
    
    % Apply perturbations separately to real and imaginary parts
    perturbedAmp = real(dipSubset.complAmpl) .* realPerturbationFactor + ...
                   1i * imag(dipSubset.complAmpl) .* imagPerturbationFactor;
    
    % Normalize perturbed dipole amplitudes
    perturbedAmpReal = real(perturbedAmp);
    perturbedAmpImag = imag(perturbedAmp);
    normalizedPerturbedAmpReal = perturbedAmpReal / max(abs(perturbedAmpReal));
    normalizedPerturbedAmpImag = perturbedAmpImag / max(abs(perturbedAmpImag));
    
    % Create an initial guess based on the perturbed and normalized amplitudes
    initialGuess = [normalizedPerturbedAmpReal; normalizedPerturbedAmpImag]';
    
    % Bounds for PSO (real and imaginary parts separately)
    realLowerBounds = min(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
    realUpperBounds = max(normalizedPerturbedAmpReal) * ones(numDipoles, 1);
    imagLowerBounds = min(normalizedPerturbedAmpImag) * ones(numDipoles, 1);
    imagUpperBounds = max(normalizedPerturbedAmpImag) * ones(numDipoles, 1);
    lB = [realLowerBounds; imagLowerBounds];
    uB = [realUpperBounds; imagUpperBounds];
    
    % Optimization function for PSO
    optimFun = @(amp) optimization.normObjectiveFunction(amp(1:numDipoles).' + ...
        1i * amp((numDipoles+1):(numDipoles*2)).', dipSubset, f0List, points, weights, fF_ref, totalPower_ref);
    
    % PSO options
    swarmSize = numDipoles * 2;  % Number of particles (twice the number of dipoles)
    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', swarmSize, ...               % Number of particles in the swarm
        'MaxIterations', 200, ...                 % Maximum number of iterations
        'InertiaRange', [0.3, 1.5], ...           % Range for inertia weight (balance exploration/exploitation)
        'SelfAdjustmentWeight', 1.1, ...          % Weight for a particle's own best experience
        'SocialAdjustmentWeight', 1.05, ...       % Weight for following the global best solution
        'FunctionTolerance', 1e-10, ...           % Stop if function value improvement is below this threshold
        'MaxStallIterations', 20, ...             % Stop if no improvement in 40 consecutive iterations
        'InitialSwarmMatrix', repmat(initialGuess, swarmSize, 1),... % Use initial guess
        'Display', 'iter' ...                      % Show progress (can be 'iter' for more detailed display)
        );
    
    % Run PSO for Amplitude Recovery
    tic;
    [optAmps_pso_vec, finalError_pso] = particleswarm(optimFun, 2 * numDipoles, lB, uB, options_pso);
    pso_time = toc;
    
    % Reconstruct the complex amplitudes
    optAmps_pso = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec(numDipoles+1:end);
    
    % Store results
    timeResults(i) = pso_time;
    finalErrorResults(i) = finalError_pso;
end

% Plot results: Time vs. Number of Dipoles
figure;
subplot(2,1,1);
plot(dipoleCounts, timeResults, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Computation Time (s)');
title('PSO Optimization Time vs. Number of Dipoles');
grid on;

% Plot results: Error vs. Number of Dipoles
subplot(2,1,2);
plot(dipoleCounts, finalErrorResults, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Final Error');
title('PSO Final Error vs. Number of Dipoles');
grid on;

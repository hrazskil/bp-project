clc; clear; close all;

% Dipole counts to test
NdipList = [100, 200, 300, 400, 500, 600, 700, 780];  % example subset counts, adapt as needed

% Load observation data (fixed)
Nleb = 170;
vfile = sprintf('tests/test_structure_2/dataHV_N_obs_points/dataHV_%d_obs_points_Vertical.tsv', Nleb);
hfile = sprintf('tests/test_structure_2/dataHV_N_obs_points/dataHV_%d_obs_points_Horizontal.tsv', Nleb);
load('tests/test_structure_2/DipoleArray.mat');  % full dipole config

dataVertical = load(vfile);
dataHorizontal = load(hfile);
dipoleFull = dip;  % full original dipole array
frequency = f0List;

% Observation points
Nphi = size(dataVertical,1);
Nbeta = size(dataHorizontal,1);
x1 = cos(dataVertical(:,1));  y1 = zeros(Nphi,1);  z1 = sin(dataVertical(:,1));
x2 = cos(dataHorizontal(:,1)); y2 = sin(dataHorizontal(:,1)); z2 = zeros(Nbeta,1);
vertical_points = [x1, y1, z1];
horizontal_points = [x2, y2, z2];

% Create input data structure
inputData.vertical.rad = dataVertical(:,2);
inputData.vertical.points = vertical_points;
inputData.vertical.weights = ones(Nphi, 1);
inputData.vertical.totalPower = sum(inputData.vertical.rad);

inputData.horizontal.rad = dataHorizontal(:,2);
inputData.horizontal.points = horizontal_points;
inputData.horizontal.weights = ones(Nbeta, 1);
inputData.horizontal.totalPower = sum(inputData.horizontal.rad);

inputData.freq = frequency;

% Physical constants
construct = utilities.constants.giveConstants();
omega = 2 * pi * frequency;
k = omega / construct.c0;
rFar = 1e6 / k;

[points, weights, ~] = utilities.getLebedevSphere(Nleb);
rObserved = points * rFar;

% Reference full far-field pattern
fF_ref_full = fieldEvaluation.farFieldM2(rObserved, dipoleFull, frequency);
totalPower_ref = sum(sum(fF_ref_full .* conj(fF_ref_full), 2) .* weights) / (2 * construct.Z0);
dipoleFull.complAmpl = dipoleFull.complAmpl / sqrt(totalPower_ref);

% Results struct
results(length(NdipList)) = struct('Ndip', [], 'optAmps_fmincon', [], 'perturbationFactor', [], 'selectedIdx', []);

for i = 1:length(NdipList)
    numDipoles = NdipList(i);

    % Subset dipole model using evenly spread indices
    selectedIdx = round(linspace(1, numel(dipoleFull.complAmpl), numDipoles));
    dipSubset = struct();
    dipSubset.complAmpl = dipoleFull.complAmpl(selectedIdx);
    dipSubset.pos = dipoleFull.pos(selectedIdx, :);
    dipSubset.dir = dipoleFull.dir(selectedIdx, :);

    % Perturb subset
    perturbationFactor = 0.1 * abs(dipSubset.complAmpl) .* ...
                         (randn(numDipoles,1) + 1i * randn(numDipoles,1));
    dipPert = dipSubset;
    dipPert.complAmpl = dipSubset.complAmpl + perturbationFactor;

    % Normalize perturbed dipole
    fF_pert = fieldEvaluation.farFieldM2(rObserved, dipPert, frequency);
    totalPower_pert = sum(sum(fF_pert .* conj(fF_pert), 2) .* weights) / (2 * construct.Z0);
    dipPert.complAmpl = dipPert.complAmpl / sqrt(totalPower_pert);

    % Optimization bounds
    realPert = real(dipPert.complAmpl);
    imagPert = imag(dipPert.complAmpl);
    lB = [min(realPert)*ones(numDipoles,1); min(imagPert)*ones(numDipoles,1)];
    uB = [max(realPert)*ones(numDipoles,1); max(imagPert)*ones(numDipoles,1)];

    % PSO setup
    initialGuess = [real(dipPert.complAmpl); imag(dipPert.complAmpl)]';
    swarmSize = ceil(780/4);
    initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);

    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', swarmSize, ...
        'MaxIterations', 200, ...
        'InertiaRange', [0.3, 1.5], ...
        'SelfAdjustmentWeight', 1.1, ...
        'SocialAdjustmentWeight', 1.05, ...
        'FunctionTolerance', 1e-6, ...
        'MaxStallIterations', 20, ...
        'InitialSwarmMatrix', initialSwarmMatrix, ...
        'Display', 'iter');

    optimFun = @(amp) optimization.optimFunX(amp, dipSubset, inputData, numDipoles);
    [optAmps_pso_vec, ~] = particleswarm(optimFun, 2*numDipoles, lB, uB, options_pso);
    optAmps_pso = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec(numDipoles+1:end);

    % fmincon refinement
    initialGuess = [real(optAmps_pso); imag(optAmps_pso)]';
    options_fmincon = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 50000, ...
        'OptimalityTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'Display', 'none');

    [optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun, ...
        initialGuess, [], [], [], [], lB, uB, [], options_fmincon);
    optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                      1i * optAmps_fmincon_vec(numDipoles+1:end).';

    % Store results
    results(i).Ndip = numDipoles;
    results(i).optAmps_fmincon = optAmps_fmincon;
    results(i).perturbationFactor = perturbationFactor;
    results(i).selectedIdx = selectedIdx;

    fprintf('Dipoles = %d, Final Error: %.3e\n', numDipoles, finalError_fmincon);
end

% Save result
save('results/test_structure_2_dipole_scaling_results2.mat', 'results');
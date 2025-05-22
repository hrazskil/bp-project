clc; clear; close all;

% Set of observation point counts
NobsList = [50, 170, 302, 590, 770];

% Load dipole geometry and reference frequency
load('tests/test_structure_2/DipoleArray.mat');  % loads 'dip' and 'f0List'
dipoleRef = dip;
frequency = f0List;
numDipoles = numel(dip.complAmpl);

% Create struct arrays to store results
results(length(NobsList)) = struct('Nobs', [], 'optAmps_fmincon', [], 'perturbationFactor', []);

for i = 1:length(NobsList)
    Nobs = NobsList(i);

    % Construct file paths
    vfile = sprintf('tests/test_structure_2/dataHV_N_obs_points/dataHV_%d_obs_points_Vertical.tsv', Nobs);
    hfile = sprintf('tests/test_structure_2/dataHV_N_obs_points/dataHV_%d_obs_points_Horizontal.tsv', Nobs);

    % Load far-field data
    dataVertical = load(vfile);
    dataHorizontal = load(hfile);

    % Construct observation points
    Nphi = size(dataVertical,1);
    Nbeta = size(dataHorizontal,1);

    x1 = cos(dataVertical(:,1));  y1 = zeros(Nphi,1);  z1 = sin(dataVertical(:,1));
    x2 = cos(dataHorizontal(:,1)); y2 = sin(dataHorizontal(:,1)); z2 = zeros(Nbeta,1);

    vertical_points = [x1, y1, z1];
    horizontal_points = [x2, y2, z2];

    % Construct inputData structure
    inputData.horizontal.rad        = dataHorizontal(:,2);
    inputData.horizontal.points     = horizontal_points;
    inputData.horizontal.weights    = ones(Nbeta, 1);
    inputData.horizontal.totalPower = sum(inputData.horizontal.rad);

    inputData.vertical.rad          = dataVertical(:,2);
    inputData.vertical.points       = vertical_points;
    inputData.vertical.weights      = ones(Nphi, 1);
    inputData.vertical.totalPower   = sum(inputData.vertical.rad);

    inputData.freq = frequency;

    % Perturbation
    perturbationFactor = 0.1 * abs(dipoleRef.complAmpl) .* ...
                         (randn(numDipoles,1) + 1i*randn(numDipoles,1));
    dipolePerturbed = dipoleRef;
    dipolePerturbed.complAmpl = dipoleRef.complAmpl + perturbationFactor;

    % Normalization constants
    construct = utilities.constants.giveConstants();
    omega = 2 * pi * frequency;
    k = omega / construct.c0;
    rFar = 1e6 / k;

    Nleb = 302;
    [points, weights, ~] = utilities.getLebedevSphere(Nleb);
    rObserved = points * rFar;

    fF_ref = fieldEvaluation.farFieldM2(rObserved, dipoleRef, frequency);
    fF_Pert = fieldEvaluation.farFieldM2(rObserved, dipolePerturbed, frequency);

    totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);
    totalPower_pert = sum(sum(fF_Pert .* conj(fF_Pert), 2) .* weights) / (2 * construct.Z0);

    dipoleRef.complAmpl = dipoleRef.complAmpl / sqrt(totalPower_ref);
    dipolePerturbed.complAmpl = dipolePerturbed.complAmpl / sqrt(totalPower_pert);

    % Define bounds
    realPert = real(dipolePerturbed.complAmpl);
    imagPert = imag(dipolePerturbed.complAmpl);
    lB = [min(realPert)*ones(numDipoles,1); min(imagPert)*ones(numDipoles,1)];
    uB = [max(realPert)*ones(numDipoles,1); max(imagPert)*ones(numDipoles,1)];

    % Initial guess
    initialGuess = [real(dipolePerturbed.complAmpl); imag(dipolePerturbed.complAmpl)]';

    % PSO Optimization
    swarmSize = numDipoles / 10;
    initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);

    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', swarmSize, ...
        'MaxIterations', 200, ...
        'InertiaRange', [0.3, 1.5], ...
        'SelfAdjustmentWeight', 1.1, ...
        'SocialAdjustmentWeight', 1.05, ...
        'FunctionTolerance', 1e-5, ...
        'MaxStallIterations', 20, ...
        'InitialSwarmMatrix', initialSwarmMatrix, ...
        'Display', 'iter');
    warning('off', 'optimlib:fwdFinDiffTooClose');
    optimFun = @(amp) optimization.optimFunX(amp, dip, inputData, numDipoles);
    [optAmps_pso_vec, ~] = particleswarm(optimFun, 2*numDipoles, lB, uB, options_pso);
    optAmps_pso = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec(numDipoles+1:end);
    warning('on', 'optimlib:fwdFinDiffTooClose');
    % fmincon Refinement
    initialGuess = [real(optAmps_pso); imag(optAmps_pso)]';
    options_fmincon = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'MaxIterations', 50, ...
        'MaxFunctionEvaluations', 10000, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'Display', 'none');

    [optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun, ...
        initialGuess, [], [], [], [], lB, uB, [], options_fmincon);

    optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                      1i * optAmps_fmincon_vec(numDipoles+1:end).';

    % Save results
    results(i).Nobs = Nobs;
    results(i).optAmps_fmincon = optAmps_fmincon;
    results(i).perturbationFactor = perturbationFactor;

    fprintf('N = %d, Final Error: %.3e\n', Nobs, finalError_fmincon);
end



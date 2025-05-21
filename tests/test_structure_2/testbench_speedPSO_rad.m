%% --- PSO Optimization with normObjectiveFunction_rad ---
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


% Set test dipole subset sizes
dipoleCounts = round(linspace(1, numel(dip.complAmpl)/4, 4));

% Prepare outputs
timeResults = zeros(length(dipoleCounts), 1);
finalErrorResults = zeros(length(dipoleCounts), 1);

for i = 1:length(dipoleCounts)
    numDipoles = dipoleCounts(i);
    dipSubset.complAmpl = dip.complAmpl(1:numDipoles);
    dipSubset.pos = dip.pos(1:numDipoles,:);
    dipSubset.dir = dip.dir(1:numDipoles,:);

    % Perturbation
    realPerturb = 1 + 0.01 * randn(numDipoles,1);
    imagPerturb = 1 + 0.01 * randn(numDipoles,1);
    pertAmp = real(dipSubset.complAmpl).*realPerturb + 1i*imag(dipSubset.complAmpl).*imagPerturb;

    normReal = real(pertAmp) / max(abs(real(pertAmp)));
    normImag = imag(pertAmp) / max(abs(imag(pertAmp)));

    initialGuess = [normReal; normImag]';

    % Bounds
    lB = [min(normReal)*ones(numDipoles,1); min(normImag)*ones(numDipoles,1)];
    uB = [max(normReal)*ones(numDipoles,1); max(normImag)*ones(numDipoles,1)];

    % Define objective function
    optimFun = @(amp) optimization.normObjectiveFunction_rad(...
        updateDipoleAmpl(dipSubset, amp(1:numDipoles).' + 1i*amp(numDipoles+1:end).'), ...
        inputData);

    % PSO settings
    swarmSize = numDipoles * 2;
    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', swarmSize, ...
        'MaxIterations', 100, ...
        'Display', 'iter', ...
        'FunctionTolerance', 1e-5, ...
        'MaxStallIterations', 30, ...
        'InitialSwarmMatrix', repmat(initialGuess, swarmSize, 1));

    % Run PSO
    tic;
    [optVec, err] = particleswarm(optimFun, swarmSize, lB, uB, options_pso);
    timeResults(i) = toc;
    finalErrorResults(i) = err;
end

% Plot: Time vs. Dipole Count
figure;
subplot(2,1,1);
plot(dipoleCounts, timeResults, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles'); ylabel('Computation Time (s)');
title('PSO Time vs. Number of Dipoles'); grid on;

% Plot: Error vs. Dipole Count
subplot(2,1,2);
plot(dipoleCounts, finalErrorResults, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles'); ylabel('Final Error');
title('PSO Final Error vs. Number of Dipoles'); grid on;

% --- Helper function ---
function dipOut = updateDipoleAmpl(dipIn, newAmp)
    dipOut = dipIn;
    dipOut.complAmpl = newAmp;
end

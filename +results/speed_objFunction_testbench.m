%% --- Benchmarking normObjectiveFunction ---
clc; clear; close all;

% Load dipole variables from file 
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\DipoleArray.mat');

% Extract dipole information
numDipoles = numel(dip.complAmpl);

% Define available Lebedev quadrature degrees
lebedevDegrees = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810];
numPointsRange = lebedevDegrees; % Use only available Lebedev degrees

% Define physical constants
construct = utilities.constants.giveConstants();
omega = 2 * pi * f0List;  % Angular frequency
k = omega / construct.c0;    % Wavenumber
rFar = 1e6 / k;              % Large observation distance

% Generate test data for all Lebedev quadratures
[pointsCell, weightsCell] = arrayfun(@(numPoints) utilities.getLebedevSphere(numPoints), numPointsRange, 'UniformOutput', false);
rObservedCell = cellfun(@(points) points * rFar, pointsCell, 'UniformOutput', false);

% Compute reference far-field patterns and total power
fF_refCell = cellfun(@(rObs) fieldEvaluation.farField(rObs, dip, f0List), rObservedCell, 'UniformOutput', false);
totalPower_refCell = cellfun(@(fF_ref, weights) sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0), fF_refCell, weightsCell);

% Time measurement for all cases
timeResults = cellfun(@(points, weights, fF_ref, totalPower_ref) ...
    timeit(@() optimization.normObjectiveFunction(dip.complAmpl, dip, f0List, points, weights, fF_ref, totalPower_ref)), ...
    pointsCell, weightsCell, fF_refCell, num2cell(totalPower_refCell));

% Plot results
figure;
semilogx(numPointsRange, timeResults, '-o', 'LineWidth', 2);
xlabel('Number of Observation Points (Lebedev Degrees)');
ylabel('Computation Time (s)');
title('Benchmark of normObjectiveFunction with Varying degrees');
grid on;

%% --- Benchmarking Different Numbers of Dipoles (Fixed Lebedev 302) ---
dipoleCounts = round(linspace(10, numDipoles, 200)); % Different dipole counts
degree = 302; % Fixed Lebedev quadrature degree

% Generate Lebedev quadrature points and weights
[points, weights, ~] = utilities.getLebedevSphere(degree);
rObserved = points * rFar;   % Scale points to observation distance

% Compute reference far-field pattern for fixed Lebedev 302
fF_ref = fieldEvaluation.farField(rObserved, dip, f0List);
totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);

% Generate dipole subsets using arrayfun for vectorization
complAmplSubsets = arrayfun(@(n) dip.complAmpl(1:n), dipoleCounts, 'UniformOutput', false);
positionsSubsets = arrayfun(@(n) dip.pos(1:n, :), dipoleCounts, 'UniformOutput', false);
directionsSubsets = arrayfun(@(n) dip.dir(1:n, :), dipoleCounts, 'UniformOutput', false);

% Create dipole subset structs using cellfun
dipSubsets = cellfun(@(complAmpl, pos, dir) struct('complAmpl', complAmpl, ...
                                                   'pos', pos, ...
                                                   'dir', dir), ...
    complAmplSubsets, positionsSubsets, directionsSubsets, 'UniformOutput', false);

% Time measurement for all dipole subsets using arrayfun for vectorization
timeResultsDipoles = arrayfun(@(i) timeit(@() optimization.normObjectiveFunction(dipSubsets{i}.complAmpl, dipSubsets{i}, f0List, points, weights, fF_ref, totalPower_ref)), 1:length(dipSubsets));

% Plot results
figure;
plot(dipoleCounts, timeResultsDipoles, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Computation Time (s)');
title('normObjectiveFunction with Varying Dipoles (Lebedev 302)');
grid on;
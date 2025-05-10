%% --- Benchmarking normObjectiveFunction ---
clc; clear; close all;

% Load dipole variables from file 
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat')

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

% --- Generate test data for all Lebedev quadratures ---
[pointsCell, weightsCell] = arrayfun(@(numPoints) utilities.getLebedevSphere(numPoints), numPointsRange, 'UniformOutput', false);
rObservedCell = cellfun(@(points) points * rFar, pointsCell, 'UniformOutput', false);

% --- Initialize output arrays ---
numRepeats = 20;
timeResults = zeros(size(numPointsRange));
fF_refCell = cell(size(numPointsRange));
totalPower_refCell = cell(size(numPointsRange));

% --- Compute far-fields and benchmark with repeated timing ---
for i = 1:length(numPointsRange)
    points = pointsCell{i};
    weights = weightsCell{i};
    rObserved = rObservedCell{i};

    % Compute reference far-field and total power
    fF_ref = fieldEvaluation.farFieldM2(rObserved, dip, f0List);
    totalPower = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);

    % Store for potential reuse
    fF_refCell{i} = fF_ref;
    totalPower_refCell{i} = totalPower;

    % Benchmark repeated evaluations
    times = zeros(numRepeats, 1);
    for r = 1:numRepeats
        times(r) = timeit(@() optimization.normObjectiveFunction(dip.complAmpl, dip, f0List, points, weights, fF_ref, totalPower));
    end
    timeResults(i) = mean(times);
end


% Plot results
figure;
semilogx(numPointsRange, timeResults, '-o', 'LineWidth', 2);
xlabel('Number of Observation Points (Lebedev Degrees)');
ylabel('Computation Time (s)');
title('Benchmark of normObjectiveFunction with Varying degrees');
grid on;

%% --- Recompute timings for specific Lebedev degrees ---
degreesToRecompute = [5810];  % <- Replace with any degrees you'd like to redo
numRepeats = 3;                   % Increase for better averaging

% Loop over each requested degree
for deg = degreesToRecompute
    idx = find(lebedevDegrees == deg);
    if isempty(idx)
        warning('Degree %d not found in lebedevDegrees array.', deg);
        continue;
    end

    % Retrieve precomputed components
    points = pointsCell{idx};
    weights = weightsCell{idx};
    rObserved = rObservedCell{idx};
    fF_ref = fF_refCell{idx};
    totalPower = totalPower_refCell{idx};

    % Time the objective function
    times = zeros(numRepeats, 1);
    for r = 1:numRepeats
        times(r) = timeit(@() optimization.normObjectiveFunction(...
            dip.complAmpl, dip, f0List, points, weights, fF_ref, totalPower));
    end
    timeResults(idx) = mean(times);
    fprintf('Recomputed time for degree %d (index %d): %.4f s\n', ...
        deg, idx, timeResults(idx));
end

%% --- Benchmarking Different Numbers of Dipoles (Fixed Lebedev 302) ---
dipoleCounts = round(linspace(10, numDipoles, 40)); % Different dipole counts
degree = 770; % Fixed Lebedev quadrature degree

% Generate Lebedev quadrature points and weights
[points, weights, ~] = utilities.getLebedevSphere(degree);
rObserved = points * rFar;   % Scale points to observation distance

% Compute reference far-field pattern for fixed Lebedev 302
fF_ref = fieldEvaluation.farFieldM2(rObserved, dip, f0List);
totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);

% Initialize result array
timeResultsDipoles = zeros(size(dipoleCounts));
numRepeats = 100;

% Loop over dipole counts
for i = 1:length(dipoleCounts)
    n = dipoleCounts(i);
    
    % Create subset dipole structure
    subDip.complAmpl = dip.complAmpl(1:n);
    subDip.pos = dip.pos(1:n, :);
    subDip.dir = dip.dir(1:n, :);
    
    % Repeat timing
    tvals = zeros(numRepeats, 1);
    for r = 1:numRepeats
        tvals(r) = timeit(@() optimization.normObjectiveFunction(subDip.complAmpl, subDip, f0List, points, weights, fF_ref, totalPower_ref));
    end
    timeResultsDipoles(i) = mean(tvals);
    dipoleCounts(i)
end

% Plot results
figure;
plot(dipoleCounts, timeResultsDipoles, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Computation Time (s)');
title('normObjectiveFunction with Varying Dipoles (Lebedev 302)');
grid on;


%% --- Recompute timings for specific dipole counts ---
dipoleCountsToRecompute = [681];  % <- Replace with dipole counts to retest
numRepeats = 5;                       % Averaging repetitions

for dipVal = dipoleCountsToRecompute
    idx = find(dipoleCounts == dipVal);
    if isempty(idx)
        warning('Dipole count %d not found in dipoleCounts array.', dipVal);
        continue;
    end

    % Extract subset of dipoles
    n = dipoleCounts(idx);
    subDip.complAmpl = dip.complAmpl(1:n);
    subDip.pos = dip.pos(1:n, :);
    subDip.dir = dip.dir(1:n, :);

    % Recompute timing
    tvals = zeros(numRepeats, 1);
    for r = 1:numRepeats
        tvals(r) = timeit(@() optimization.normObjectiveFunction( ...
            subDip.complAmpl, subDip, f0List, points, weights, fF_ref, totalPower_ref));
    end
    timeResultsDipoles(idx) = mean(tvals);
    fprintf('Recomputed time for %d dipoles (index %d): %.4f s\n', ...
        n, idx, timeResultsDipoles(idx));
end

% Optional: Refresh plot
figure;
plot(dipoleCounts, timeResultsDipoles, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Computation Time (s)');
title('Updated normObjectiveFunction Time vs. Number of Dipoles');
grid on;

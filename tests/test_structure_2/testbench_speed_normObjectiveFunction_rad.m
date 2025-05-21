%% --- Benchmarking normObjectiveFunction_rad ---
clc; clear; close all;

%% Load test dipole configuration and reference far-field data
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat');

%% Save generated far-field test data for test_structure_2
N = 360; % number of Points
angles = linspace(0, 2*pi*(360-1)/360, N)';

% Horizontal plane: xy
horizontal_points = [cos(angles), sin(angles), zeros(N,1)];
E_horiz = fieldEvaluation.farFieldM2(horizontal_points, dip, f0List);
I_horiz = fieldEvaluation.powerDensityFar(E_horiz);

% Vertical plane: xz
vertical_points = [cos(angles), zeros(N,1), sin(angles)];
E_vert = fieldEvaluation.farFieldM2(vertical_points, dip, f0List);
I_vert = fieldEvaluation.powerDensityFar(E_vert);

% Save horizontal
dataH = [angles, I_horiz];
fidH = fopen('C:\\Users\\kilia\\Plocha\\gitHub\\bp-project\\tests\\test_structure_2\\dataHorizontal.tsv', 'w');
fprintf(fidH, '%.16e\t%.16e\n', dataH.');
fclose(fidH);

% Save vertical
dataV = [angles, I_vert];
fidV = fopen('C:\\Users\\kilia\\Plocha\\gitHub\\bp-project\\tests\\test_structure_2\\dataVertical.tsv', 'w');
fprintf(fidV, '%.16e\t%.16e\n', dataV.');
fclose(fidV);

disp('Far-field .tsv files saved successfully for test_structure_2.');

%% Construct inputData struct from vertical and horizontal scans
Nphi  = size(dataV,1);
Nbeta = size(dataH,1);

x1 = cos(dataV(:,1)); y1 = zeros(Nphi,1); z1 = sin(dataV(:,1));
x2 = cos(dataH(:,1)); y2 = sin(dataH(:,1)); z2 = zeros(Nbeta,1);

inputData.vertical.points       = [x1 y1 z1];
inputData.vertical.rad          = dataV(:,2);
inputData.vertical.weights      = ones(Nphi, 1);
inputData.vertical.totalPower   = sum(inputData.vertical.rad .* inputData.vertical.weights);

inputData.horizontal.points     = [x2 y2 z2];
inputData.horizontal.rad        = dataH(:,2);
inputData.horizontal.weights    = ones(Nbeta, 1);
inputData.horizontal.totalPower = sum(inputData.horizontal.rad .* inputData.horizontal.weights);

inputData.freq = f0List;  % Assign frequency list

%% Evaluate reference error with normObjectiveFunction_rad
tic
error = optimization.normObjectiveFunction_rad(dip, inputData);
toc
disp(['Initial Error (Normalized): ', num2str(error)]);

%% Benchmark runtime vs number of dipoles
dipoleCounts = round(linspace(10, size(dip.complAmpl,1), 40));
timeResults = zeros(size(dipoleCounts));

for i = 1:length(dipoleCounts)
    n = dipoleCounts(i);
    subDip = dip;
    subDip.complAmpl = dip.complAmpl(1:n);
    subDip.pos = dip.pos(1:n,:);
    subDip.dir = dip.dir(1:n,:);

    trials = 5;
    times = zeros(trials,1);
    for t = 1:trials
        times(t) = timeit(@() optimization.normObjectiveFunction_rad(subDip, inputData));
    end
    timeResults(i) = mean(times);
    dipoleCounts(i)
end

%% Plot benchmark results
figure;
plot(dipoleCounts, timeResults, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Computation Time (s)');
title('Runtime of normObjectiveFunction\_rad vs Dipole Count');
grid on;

%% --- Recompute timings for specific dipole counts (normObjectiveFunction_rad) ---
dipoleCountsToRecompute = [69,405,760];  % <- Specify desired dipole counts to recompute
numRepeats = 100;                       % Number of repetitions for averaging

for dipVal = dipoleCountsToRecompute
    idx = find(dipoleCounts == dipVal);
    if isempty(idx)
        warning('Dipole count %d not found in dipoleCounts array.', dipVal);
        continue;
    end

    % Extract subset dipole configuration
    n = dipoleCounts(idx);
    subDip.complAmpl = dip.complAmpl(1:n);
    subDip.pos = dip.pos(1:n, :);
    subDip.dir = dip.dir(1:n, :);

    % Recalculate execution time
    tvals = zeros(numRepeats, 1);
    for r = 1:numRepeats
        tvals(r) = timeit(@() optimization.normObjectiveFunction_rad(subDip, inputData));
    end

    timeResults(idx) = mean(tvals);
    fprintf('Recomputed time for %d dipoles (index %d): %.4f s\n', ...
        n, idx, timeResults(idx));
end

% Optional: Refresh plot
figure;
plot(dipoleCounts, timeResults, '-o', 'LineWidth', 2);
xlabel('Number of Dipoles');
ylabel('Computation Time (s)');
title('Updated Runtime of normObjectiveFunction\_rad vs Dipole Count');
grid on;


%% Benchmark runtime vs number of observation points
pointCounts = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]./2;
timePerPoints = zeros(size(pointCounts));

for i = 1:length(pointCounts)
    N = pointCounts(i);
    angles = linspace(0, 2*pi*(N-1)/N, N)';

    vertical_points = [cos(angles), zeros(N,1), sin(angles)];
    horizontal_points = [cos(angles), sin(angles), zeros(N,1)];
    size(vertical_points,1)
    E_vert = fieldEvaluation.farFieldM2(vertical_points, dip, f0List);
    I_vert = fieldEvaluation.powerDensityFar(E_vert);
    E_horiz = fieldEvaluation.farFieldM2(horizontal_points, dip, f0List);
    I_horiz = fieldEvaluation.powerDensityFar(E_horiz);

    input.vertical.points       = vertical_points;
    input.vertical.rad          = I_vert;
    input.vertical.weights      = ones(N,1);
    input.vertical.totalPower   = sum(I_vert);

    input.horizontal.points     = horizontal_points;
    input.horizontal.rad        = I_horiz;
    input.horizontal.weights    = ones(N,1);
    input.horizontal.totalPower = sum(I_horiz);

    input.freq = f0List;

    trials = 20;
    tval = zeros(trials,1);
    for t = 1:trials
        tval(t) = timeit(@() optimization.normObjectiveFunction_rad(dip, input));
    end
    timePerPoints(i) = mean(tval);
end

%% Plot runtime vs number of observation points
figure;
plot(pointCounts*2, timePerPoints, '-s', 'LineWidth', 2);
xlabel('Number of Observation Points (both planes sumed)');
ylabel('Computation Time (s)');
title('Runtime of normObjectiveFunction\_rad vs Number of Points');
grid on;

%% --- Recompute timings for specific numbers of observation points ---
pointCountsToRecompute = [5810]/2;  % <- Half the total point count (per plane) (insert total number)
numRepeats = 3;

for ptVal = pointCountsToRecompute
    idx = find(pointCounts == ptVal);
    if isempty(idx)
        warning('Point count %d not found in pointCounts array.', ptVal);
        continue;
    end

    N = pointCounts(idx);
    angles = linspace(0, 2*pi*(N-1)/N, N)';

    vertical_points = [cos(angles), zeros(N,1), sin(angles)];
    horizontal_points = [cos(angles), sin(angles), zeros(N,1)];

    E_vert = fieldEvaluation.farFieldM2(vertical_points, dip, f0List);
    I_vert = fieldEvaluation.powerDensityFar(E_vert);

    E_horiz = fieldEvaluation.farFieldM2(horizontal_points, dip, f0List);
    I_horiz = fieldEvaluation.powerDensityFar(E_horiz);

    input.vertical.points       = vertical_points;
    input.vertical.rad          = I_vert;
    input.vertical.weights      = ones(N,1);
    input.vertical.totalPower   = sum(I_vert);

    input.horizontal.points     = horizontal_points;
    input.horizontal.rad        = I_horiz;
    input.horizontal.weights    = ones(N,1);
    input.horizontal.totalPower = sum(I_horiz);

    input.freq = f0List;

    % Re-timing
    tvals = zeros(numRepeats, 1);
    for r = 1:numRepeats
        tvals(r) = timeit(@() optimization.normObjectiveFunction_rad(dip, input));
    end
    timePerPoints(idx) = mean(tvals);
    fprintf('Recomputed time for %d total points (index %d): %.4f s\n', ...
        2*N, idx, timePerPoints(idx));
end

% Optional: Update plot
figure;
plot(pointCounts*2, timePerPoints, '-s', 'LineWidth', 2);
xlabel('Number of Observation Points (both planes summed)');
ylabel('Computation Time (s)');
title('Updated Runtime of normObjectiveFunction\_rad vs Number of Points');
grid on;

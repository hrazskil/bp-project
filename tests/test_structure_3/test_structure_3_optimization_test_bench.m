%% --- Initialization for BTS data without DipoleArray.mat ---
clc; clear; close all;

%% Load saved far-field data
load('tests/test_structure_3/BTSdataVertical.tsv');   % vertical far-field data
load('tests/test_structure_3/BTSdataHorizontal.tsv'); % horizontal far-field data

%% Number of observation points
Nbeta = size(BTSdataHorizontal(:,2),1);  % number of azimuthal observations
Nphi = size(BTSdataVertical(:,2),1);     % number of elevation observations

% Reconstruct Cartesian observation points from angles
x1 = cos(BTSdataVertical(:,1));
y1 = 0*x1;
z1 = sin(BTSdataVertical(:,1));

x2 = cos(BTSdataHorizontal(:,1));
y2 = sin(BTSdataHorizontal(:,1));
z2 = 0*x2;

vertical_points = [x1, y1, z1];
horizontal_points = [x2, y2, z2];

% Create input data structure for optimization
inputData.horizontal.rad = BTSdataHorizontal(:,2);
inputData.horizontal.points = horizontal_points;
inputData.horizontal.weights = ones(Nbeta, 1);
inputData.horizontal.totalPower = sum(inputData.horizontal.rad);

inputData.vertical.rad = BTSdataVertical(:,2);
inputData.vertical.points = vertical_points;
inputData.vertical.weights = ones(Nphi, 1);
inputData.vertical.totalPower = sum(inputData.vertical.rad);

inputData.freq = 785000000.000000;

%% Dipole grid generation
% === Parameters for Crossed Dipole "Pluses" ===
numPluses = 4;                    % Total number of pluses
dipoleLength = 1;              % Length of each dipole [m]
numDipolesPerArm = 20;            % Dipoles per arm (X and Z)
offsetStep = 1.5;                % Spacing between pluses along Z
rotationX = pi/2;                 % Dipoles along Z (vertical arm)
rotationY = 0;                    % Dipoles along X (horizontal arm)

amplitudes = ones(numDipolesPerArm, 1);  % Uniform excitation

% === Initialize dipole structure ===
dip.pos = [];
dip.dir = [];
dip.complAmpl = [];

% === Generate Crossed Pluses ===
for i = 1:numPluses
    % Centered around Z = 0
    zOffset = (i - (numPluses + 1)/2) * offsetStep;
    origin = [0, 0, zOffset];  % All at x=0, y=0

    % Horizontal dipole line (X-direction)
    dipX = geometry.halfwaveDipoleArray(dipoleLength, numDipolesPerArm, rotationY, origin, amplitudes);

    % Vertical dipole line (Z-direction)
    dipZ = geometry.halfwaveDipoleArray(dipoleLength, numDipolesPerArm, rotationX, origin, amplitudes);

    % Append both arms to full dipole array
    dip.pos        = [dip.pos; dipX.pos; dipZ.pos];
    dip.dir        = [dip.dir; dipX.dir; dipZ.dir];
    dip.complAmpl  = [dip.complAmpl; dipX.complAmpl; dipZ.complAmpl];
end

% === Grid Parameters (YZ-plane) ===
gridSizeY = 5;                   % Number of points along Y
gridSizeZ = 5;                   % Number of points along Z
gridWidthY = 1;                % Total physical width along Y [m]
gridHeightZ = 1;               % Total physical height along Z [m]
gridOffsetX = -0.38;             % X-offset behind each plus
rotationGrid = 0;                % Dipoles oriented along X

% Compute spacing from desired total width and height
gridSpacingY = gridWidthY / (gridSizeY - 1);
gridSpacingZ = gridHeightZ / (gridSizeZ - 1);

% Amplitudes for each grid dipole
gridDipolesPerPlane = gridSizeY * gridSizeZ;
amplitudesGrid = ones(gridDipolesPerPlane, 1);

% Optional: compute height of plus to align grid
plusHeight = (numDipolesPerArm - 1) * dipoleLength;

for i = 1:numPluses
    % Z position of the plus
    zOffset = (i - (numPluses + 1)/2) * offsetStep;

    % Centered in Y and Z to match the plus geometry
    [Yg, Zg] = meshgrid( ...
        linspace(-gridWidthY/2, gridWidthY/2, gridSizeY), ...
        linspace(zOffset - gridHeightZ/2, zOffset + gridHeightZ/2, gridSizeZ) ...
    );
    Xg = ones(size(Yg)) * gridOffsetX;
    gridPos = [Xg(:), Yg(:), Zg(:)];

    % Generate dipoles
    dipGrid = geometry.halfwaveDipoleArray(dipoleLength, gridDipolesPerPlane, ...
                 rotationGrid, [0, 0, 0], amplitudesGrid);
    dipGrid.pos = gridPos;

    % Append to dipole structure
    dip.pos        = [dip.pos; dipGrid.pos];
    dip.dir        = [dip.dir; dipGrid.dir];
    dip.complAmpl  = [dip.complAmpl; dipGrid.complAmpl];
end

% === Optional Visualization ===
figure;
quiver3(dip.pos(:,1), dip.pos(:,2), dip.pos(:,3), ...
        dip.dir(:,1), dip.dir(:,2), dip.dir(:,3), ...
        0.01, 'LineWidth', 1.2);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Crossed Half-Wave Dipole Pairs (Stacked)');
grid on;

% === Save dipole structure ===
save('dipoleStructure.mat', 'dip');
disp('Dipole structure saved to halfwaveDipole.mat.');

% Extract numDipoles
numDipoles = numel(dip.complAmpl);

% Compute the maximum magnitude across all dipole amplitudes
maxAmp = max(abs(dip.complAmpl));

% Normalize the complex amplitudes
dipoleRef = dip;
dipoleRef.complAmpl = dip.complAmpl / maxAmp;

%% --- 2 Perturbation of Initial Amplitudes ---
% Generate complex perturbation factors with small variations
perturbationFactor = 1 + 0.1 * (randn(numDipoles, 1));

% Apply perturbations to the normalized complex amplitudes
perturbedAmp = dipoleRef.complAmpl .* perturbationFactor;

% Normalize the perturbed amplitudes
maxPerturbedAmp = max(abs(perturbedAmp));
dipolePerturbedRef = dip;
dipolePerturbedRef.complAmpl = perturbedAmp / maxPerturbedAmp;

% Use perturbed amplitudes for error evaluation
dip.complAmpl = dipolePerturbedRef.complAmpl;
tic
error = optimization.normObjectiveFunction_rad(dip, inputData);
toc
disp(['Initial Test Error perturbed dip: ', num2str(error)]);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% --- 3. Optimization Using PSO ---

% --- Define Bounds ---
realPert = real(dipolePerturbedRef.complAmpl);
imagPert = imag(dipolePerturbedRef.complAmpl);

ampMinReal = min(realPert);  ampMaxReal = max(realPert);
ampMinImag = min(imagPert);  ampMaxImag = max(imagPert);

lB = [ampMinReal * ones(numDipoles, 1); ampMinImag * ones(numDipoles, 1)];
uB = [ampMaxReal * ones(numDipoles, 1); ampMaxImag * ones(numDipoles, 1)];

% --- Initialize PSO Parameters ---
% Combine real and imaginary parts into a single initial guess vector
initialGuess = [real(dipolePerturbedRef.complAmpl); imag(dipolePerturbedRef.complAmpl)]';

% Create swarm matrix by replicating the initial guess (each row is a particle)
swarmSize = numDipoles/2;  % Swarm size is double the number of dipoles (real + imag)
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



simAmpReal = real(dipoleFmincon.complAmpl);
simAmpImag = imag(dipoleFmincon.complAmpl);
normalizedsimAmpReal = simAmpReal * max(abs(simAmpReal));
normalizedsimAmpImag = simAmpImag * max(abs(simAmpImag));

dipoleFmincon.complAmpl = normalizedsimAmpReal + ...
                               1i * normalizedsimAmpImag;
% Fmincon Plot
%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon, frequency, 180, 360);
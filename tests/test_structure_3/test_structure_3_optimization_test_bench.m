%% --- Initialization for BTS data without DipoleArray.mat ---
clc; clear; close all;

%% Load saved far-field data
load('tests/test_structure_3/BTSdataVertical.tsv');   % vertical far-field data
load('tests/test_structure_3/BTSdataHorizontal.tsv'); % horizontal far-field data
load('tests/test_structure_3/BTS.mat')

dipoleRef=dip;

% Determine number of observation points from data files
Nbeta = size(BTSdataHorizontal(:,2),1);  
Nphi = size(BTSdataVertical(:,2),1);    

% Reconstruct observation point vectors from saved angular values
x1 = cos(BTSdataVertical(:,1));
y1 = 0*x1;
z1 = sin(BTSdataVertical(:,1));

x2 = cos(BTSdataHorizontal(:,1));
y2 = sin(BTSdataHorizontal(:,1));
z2 = 0*x2;

vertical_points = [x1(:,1), y1(:,1), z1(:,1)];
horizontal_points = [x2(:,1), y2(:,1), z2(:,1)];

% Assemble structured input data for optimization
inputData.horizontal.rad        = BTSdataHorizontal(:,2);
inputData.horizontal.points     = horizontal_points;
inputData.horizontal.weights    = ones(Nbeta, 1);  % Uniform weights

inputData.horizontal.totalPower = ...
    sum((inputData.horizontal.rad .* inputData.horizontal.weights));

inputData.vertical.rad          = BTSdataVertical(:,2);
inputData.vertical.points       = vertical_points;
inputData.vertical.weights      = ones(Nbeta, 1);  % Uniform weights

inputData.vertical.totalPower = ...
    sum((inputData.vertical.rad .* inputData.vertical.weights));

inputData.freq = f0List;  % Frequency list used in computation


tic
    error = optimization.normObjectiveFunction_rad(dipoleRef, inputData);
    toc
    disp(['Initial Test reference to inputData: ', num2str(error)]);

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
    dip.complAmpl  = [dip.complAmpl; dipGrid.complAmpl*eps];
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

%% --- 2 Perturbation of Initial Amplitudes ---
perturbationFactor = 0.1*abs(dip.complAmpl).*(randn(numDipoles, 1)+1i*randn(numDipoles, 1));
dipolePer = dip;
dipolePer.complAmpl = dip.complAmpl + perturbationFactor;


%% Normalize dipole amplitudes
    % Define physical constants
    construct = utilities.constants.giveConstants();
    omega = 2 * pi * inputData.freq;  % Angular frequency
    k = omega / construct.c0;    % Wavenumber
    rFar = 1e6 / k;              % Large observation distance

    % degree: { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,
    %     350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,
    %     3470, 3890, 4334, 4802, 5294, 5810 };
    Nleb = 302;                 % Number of Lebedev quadrature points

    % Get Lebedev quadrature points and weights
    [points, weights, ~] = utilities.getLebedevSphere(Nleb);
    rObserved = points * rFar;   % Scale points to observation distance

    % Compute far-field patterns
    % fF_ref = fieldEvaluation.farFieldM2(rObserved, dipoleRef, inputData.freq);
    fF_Per = fieldEvaluation.farFieldM2(rObserved, dipolePer, inputData.freq);

    % Compute total radiated power for normalization
    % totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);
    totalPower_Per = sum(sum(fF_Per .* conj(fF_Per), 2) .* weights) / (2 * construct.Z0);

    % normalize
    % dipoleRef.complAmpl = dipoleRef.complAmpl / sqrt(totalPower_ref);
    dipolePer.complAmpl = dipolePer.complAmpl / sqrt(totalPower_Per);
    
    

    tic
    error = optimization.normObjectiveFunction_rad(dipolePer, inputData);
    toc
    disp(['Initial Test reference to inputData After Normalization: ', num2str(error)]);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipolePer, inputData.freq, 180, 360);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipolePer, inputData.freq, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipolePer, inputData.freq, 180, 360);

%% --- 3. Optimization Using PSO ---

% --- Define Bounds ---
realPert = real(dipolePer.complAmpl);
imagPert = imag(dipolePer.complAmpl);

ampMinReal = min(realPert);  ampMaxReal = max(realPert);
ampMinImag = min(imagPert);  ampMaxImag = max(imagPert);

lB = [ampMinReal * ones(numDipoles, 1); ampMinImag * ones(numDipoles, 1)];
uB = [ampMaxReal * ones(numDipoles, 1); ampMaxImag * ones(numDipoles, 1)];

% --- Initialize PSO Parameters ---
% Combine real and imaginary parts into a single initial guess vector
initialGuess = [real(dipolePer.complAmpl); imag(dipolePer.complAmpl)]';

% Create swarm matrix by replicating the initial guess (each row is a particle)
swarmSize = numDipoles/4;  % Swarm size is double the number of dipoles (real + imag)
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
    1i * amp(numDipoles+1:end).', dip, inputData.freq, points, weights, ...
    fF_ref, totalPower_ref);

[optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun_fmincon,...
    initialGuess, [], [], [], [], lB, uB, [], options_fmincon);

optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                  1i * optAmps_fmincon_vec(numDipoles+1:end).';
disp(['Final Error (fmincon): ', num2str(finalError_fmincon)]);


dipoleFmincon = dipoleRef;
dipoleFmincon.complAmpl = optAmps_fmincon;

% Fmincon Plot
%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360);

%% === Far-Field Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360);

%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360);
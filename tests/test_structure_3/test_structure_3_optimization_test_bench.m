%% --- Initialization for BTS data without DipoleArray.mat ---
clc; clear; close all;

%% Load saved far-field data
Fth_vert_raw   = load('tests/test_structure_3/BTSxzCutFth.tsv');   % vertical θ-cut
Fph_vert_raw   = load('tests/test_structure_3/BTSxzCutFph.tsv');   % vertical φ-cut
Fth_horiz_raw  = load('tests/test_structure_3/BTSxyCutFth.tsv');   % horizontal θ-cut
Fph_horiz_raw  = load('tests/test_structure_3/BTSxyCutFph.tsv');   % horizontal φ-cut
load('tests/test_structure_3/BTS.mat')

dipoleRef=dip;

%% === Near-Field PowerFlux: Reference ===
utilities.visualizations.PowerFluxXZ(dipoleRef, f0List, [-5.384807207388533e-01, 2.153922882955413e+00], [-5.384807207388533e-01, 5.384807207388533e-01], 201,65);

% Extract angles
theta_vert  = Fth_vert_raw(:,1);  % angle in XZ plane (elevation sweep)
phi_horiz   = Fth_horiz_raw(:,1); % angle in XY plane (azimuth sweep)

% Extract field magnitudes
Fth_vert  = Fth_vert_raw(:,2);
Fph_vert  = Fph_vert_raw(:,2);
Fth_horiz = Fth_horiz_raw(:,2);
Fph_horiz = Fph_horiz_raw(:,2);

% === Compute Total Intensity ===
intensity_vert  = Fth_vert.^2 + Fph_vert.^2;
intensity_horiz = Fth_horiz.^2 + Fph_horiz.^2;

% === Reconstruct observation directions ===
% Vertical (XZ plane): θ-sweep
x1 = sin(theta_vert); y1 = zeros(size(x1)); z1 = cos(theta_vert);
vertical_points = [x1, y1, z1];

% Horizontal (XY plane): φ-sweep
x2 = cos(phi_horiz); y2 = sin(phi_horiz); z2 = zeros(size(x2));
horizontal_points = [x2, y2, z2];

% === Assemble inputData ===
inputData.vertical.rad        = intensity_vert;
inputData.vertical.points     = vertical_points;
inputData.vertical.weights    = ones(size(intensity_vert));
inputData.vertical.totalPower = sum(intensity_vert);

inputData.horizontal.rad        = intensity_horiz;
inputData.horizontal.points     = horizontal_points;
inputData.horizontal.weights    = ones(size(intensity_horiz));
inputData.horizontal.totalPower = sum(intensity_horiz);

inputData.freq = f0List(1);  % If f0List is a vector, use first entry

tic
error = optimization.normObjectiveFunction_rad(dipoleRef, inputData);
toc
disp(['Initial Test reference to inputData: ', num2str(error)]);

%% Dipole grid generation
% Subset dipole model using evenly spread indices
numDipoles=200;
selectedIdx = round(linspace(1, numel(dipoleRef.complAmpl), numDipoles));
dipSubset = struct();
dipSubset.complAmpl = dipoleRef.complAmpl(selectedIdx);
dipSubset.pos = dipoleRef.pos(selectedIdx, :);
dipSubset.dir = dipoleRef.dir(selectedIdx, :);

indBack = dipSubset.pos(:,1) == min(dipSubset.pos(:,1));
indFront = ~indBack;
pBackMean = sum(dipSubset.complAmpl(indBack, 1))/sum(indBack);
pFrontMean = sum(dipSubset.complAmpl(indFront, 1))/sum(indFront);
dipSubset.complAmpl(indBack, 1) = pBackMean;
dipSubset.complAmpl(indFront, 1) = pFrontMean;


dipolePer = dipSubset;


% === Optional Visualization ===
% figure;
% quiver3(dipSubset.pos(:,1), dipSubset.pos(:,2), dipSubset.pos(:,3), ...
%         dipSubset.dir(:,1), dipSubset.dir(:,2), dipSubset.dir(:,3), ...
%         0.01, 'LineWidth', 1.2);
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Crossed Half-Wave Dipole Pairs (Stacked)');
% grid on;
% 
% % % === Save dipole structure ===
% % save('dipoleStructure.mat', 'dip');
% % disp('Dipole structure saved to halfwaveDipole.mat.');


error = optimization.normObjectiveFunction_rad(dipolePer, inputData);
toc
disp(['Initial Test Subset: ', num2str(error)]);

%% Normalize dipole amplitudes for Power Flux
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
% dipoleRef.complAmpl = 10*dipoleRef.complAmpl / sqrt(totalPower_ref);
dipolePer.complAmpl = 10*dipolePer.complAmpl / sqrt(totalPower_Per);



tic
error = optimization.normObjectiveFunction_rad(dipolePer, inputData);
toc
disp(['Initial Test reference to inputData After Normalization: ', num2str(error)]);



%% === Far-Field Intensity: Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef,dipolePer, f0List, 180, 360);
%% === Near-Field PowerFlux: Reference ===
utilities.visualizations.PowerFluxXZ(dipolePer, f0List, [-5.384807207388533e-01, 2.153922882955413e+00], [-5.384807207388533e-01, 5.384807207388533e-01], 201,65);

%% --- 3. Optimization Using PSO ---

dipoleRef.complAmpl = dipoleRef.complAmpl/max(abs(dipoleRef.complAmpl));
dipolePer.complAmpl = dipolePer.complAmpl/max(abs(dipolePer.complAmpl));


tic
error = optimization.normObjectiveFunction_rad(dipolePer, inputData);
toc
disp(['Initial Test reference to inputData After Normalization: ', num2str(error)]);

% --- Define Bounds ---
% realPert = real(dipolePer.complAmpl);
% imagPert = imag(dipolePer.complAmpl);
% 
% ampMinReal = min(realPert)-5;  ampMaxReal = max(realPert);
% ampMinImag = min(imagPert)-5;  ampMaxImag = max(imagPert);

% lB = [ampMinReal * ones(numDipoles, 1); ampMinImag * ones(numDipoles, 1)];
% uB = [ampMaxReal * ones(numDipoles, 1); ampMaxImag * ones(numDipoles, 1)];

lB = [-5 * ones(numDipoles, 1); -5 * ones(numDipoles, 1)];
uB = [5 * ones(numDipoles, 1); 5 * ones(numDipoles, 1)];

% --- Initialize PSO Parameters ---
% Combine real and imaginary parts into a single initial guess vector
initialGuess = [real(dipolePer.complAmpl); imag(dipolePer.complAmpl)]';

% Create swarm matrix by replicating the initial guess (each row is a particle)
    swarmSize = numDipoles*2;  % Swarm size is double the number of dipoles (real + imag)
    initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);
% swarmSize = numDipoles*2;
% % Create diverse initial swarm using uniform random sampling
% nParams = length(initialGuess);  % Should be 2*numDipoles
% 
% 
% % Uniform sampling in each dimension within [lB, uB]
% initialSwarmMatrix = rand(swarmSize, nParams) .* (uB' - lB') + lB';

    % Define PSO optimization settings
    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', swarmSize, ...                 % Number of particles
        'MaxIterations', 200, ...                   % Iteration limit
        'InertiaRange', [0.8, 1.5], ...             % Inertia control for convergence behavior
        'SelfAdjustmentWeight', 2.0, ...            % Particle's self-exploration factor
        'SocialAdjustmentWeight', 0.8, ...         % Attraction to global best solution
        'FunctionTolerance', 1e-6, ...             % Stop when improvement is below threshold
        'MaxStallIterations', 20, ...               % Stop when no progress in 40 iterations
        'InitialSwarmMatrix', initialSwarmMatrix,...% Custom initial particle positions
        'Display', 'iter' ...                       % Show iteration info in console
        );  
    
    % --- Objective Function Definition ---
    % This function evaluates the fitness of a candidate amplitude vector
    % Split real and imag parts, construct complex amplitude vector, and evaluate error
    optimFun = @(amp) optimization.optimFunX(amp, dipolePer, inputData, numDipoles);
    
    % --- Run Particle Swarm Optimization ---
    % Optimize dipole amplitudes to match far-field radiation
    [optAmps_pso_vec, finalError_pso] = particleswarm(optimFun, 2 * numDipoles, lB, uB, options_pso);
    
    
    % Combine optimized real and imaginary parts into complex amplitudes
    optAmps_pso = nan(numDipoles,1);
    optAmps_pso(:,1) = optAmps_pso_vec(1:numDipoles) + 1i * optAmps_pso_vec((numDipoles+1):(numDipoles*2));
    
    % Show final optimization error
    disp(['Final Error (PSO): ', num2str(finalError_pso)]);
    dipolePso = dipolePer;
    dipolePso.complAmpl=optAmps_pso;

%% === Near-Field PowerFlux: PSO Optimized ===
utilities.visualizations.PowerFluxXZ(dipolePso, f0List, [-0.5, 2.4], [-0.5, 0.5], 200,100);

%% --- 3. Optimization Using fmincon ---
    initialGuess = [optAmps_pso_vec(1:numDipoles);...                          % serialization of optimizations
                  optAmps_pso_vec(numDipoles+1:end)]';
    
    options_fmincon = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...             % Use Sequential Quadratic Programming (SQP) algorithm
        'MaxIterations', 50, ...           % Maximum number of iterations
        'MaxFunctionEvaluations', 1e4, ... % Maximum number of function evaluations
        'Display', 'iter');                 % Display iteration details 
    
    optimFun = @(amp) optimization.optimFunX(amp, dipSubset, inputData, numDipoles);
    
    [optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun,...
        initialGuess, [], [], [], [], lB, uB, [], options_fmincon);
    
    optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                      1i * optAmps_fmincon_vec(numDipoles+1:end).';
    disp(['Final Error (fmincon): ', num2str(finalError_fmincon)]);
    
    
    dipoleFmincon = dipolePer;
    dipoleFmincon.complAmpl = optAmps_fmincon;


%% === Far-Field Intensity Comparison: Optimized vs Reference ===
utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipoleFmincon, inputData.freq, 180, 360);


%% Normalize dipole amplitudes for Power Flux

% Compute far-field for normalization
fF_Opt = fieldEvaluation.farFieldM2(rObserved, dipoleFmincon, inputData.freq);
% Compute total radiated power for normalization
totalPower_Opt = sum(sum(fF_Opt .* conj(fF_Opt), 2) .* weights) / (2 * construct.Z0);
% normalize
dipoleFmincon.complAmpl = 10*dipoleFmincon.complAmpl / sqrt(totalPower_Opt);

%% === Near-Field PowerFlux: Result ===
utilities.visualizations.PowerFluxXZ(dipoleFmincon, f0List, [-5.384807207388533e-01, 2.153922882955413e+00], [-5.384807207388533e-01, 5.384807207388533e-01], 201,65);


% Compute far-field for normalization
fF_Opt = fieldEvaluation.farFieldM2(rObserved, dipoleRef, inputData.freq);
% Compute total radiated power for normalization
totalPower_Opt = sum(sum(fF_Opt .* conj(fF_Opt), 2) .* weights) / (2 * construct.Z0);
% normalize
dipoleRef.complAmpl = 10*dipoleRef.complAmpl / sqrt(totalPower_Opt);


utilities.visualizations.PowerFluxXAxis(dipolePer, f0List, [-5.384807207388533e-01, 2.153922882955413e+00], 201, 'dipolePerXaxis.tsv')
utilities.visualizations.PowerFluxXAxis(dipoleRef, f0List, [-5.384807207388533e-01, 2.153922882955413e+00], 201, 'dipoleRefXaxis.tsv')
utilities.visualizations.PowerFluxXAxis(dipoleFmincon, f0List, [-5.384807207388533e-01, 2.153922882955413e+00], 201, 'dipoleFminconXaxis.tsv')

type('dipolePerXaxis.tsv')
    %% --- Initialization for test_structure_2 ---
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
    dipoleRef = dip;

    % Extract numDipoles
    numDipoles = numel(dip.complAmpl);
    %% Compute initial error with true amplitudes 
    tic
    error = optimization.normObjectiveFunction_rad(dip, inputData);
    toc
    disp(['Initial Test Error: ', num2str(error)]);

    %% Create perturbed (test) amplitudes with small noise
    perturbationFactor = 0.1*abs(dipoleRef.complAmpl).*(randn(numDipoles, 1)+1i*randn(numDipoles, 1));
    dipolePerturbed = dipoleRef;
    dipolePerturbed.complAmpl = dipoleRef.complAmpl + perturbationFactor;

    % Use perturbed amplitudes for error evaluation
    dip.complAmpl = dipolePerturbed.complAmpl;
    tic
    error = optimization.normObjectiveFunction_rad(dip, inputData);
    toc
    disp(['Initial Test Error perturbed dip: ', num2str(error)]);
    %Zvláštní závislost na normalizaci

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
    fF_ref = fieldEvaluation.farFieldM2(rObserved, dipoleRef, inputData.freq);
    fF_Perturbed = fieldEvaluation.farFieldM2(rObserved, dipolePerturbed, inputData.freq);

    % Compute total radiated power for normalization
    totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);
    totalPower_Perturbed = sum(sum(fF_Perturbed .* conj(fF_Perturbed), 2) .* weights) / (2 * construct.Z0);

    % normalize
    dipoleRef.complAmpl = dipoleRef.complAmpl / sqrt(totalPower_ref);
    dipolePerturbed.complAmpl = dipolePerturbed.complAmpl / sqrt(totalPower_Perturbed);
    
    tic
    error = optimization.normObjectiveFunction_rad(dipoleRef, inputData);
    toc
    disp(['Initial Test After Normalization: ', num2str(error)]);
    % 
    % %% === Far-Field Comparison: Perturbed vs Reference ===
    % utilities.visualizations.plotFarFieldComponentComparison(dipoleRef, dipolePerturbed, inputData.freq, 180, 360);

    %% === Far-Field Intensity Comparison: Perturbed vs Reference ===
    utilities.visualizations.plotFarFieldIntensityComparison(dipoleRef, dipolePerturbed, inputData.freq, 180, 360);
    %% --- 2. Optimization Using PSO FF---
    
    % --- Define Bounds ---
    realPert = real(dipolePerturbed.complAmpl);
    imagPert = imag(dipolePerturbed.complAmpl);

    ampMinReal = min(realPert);  ampMaxReal = max(realPert);
    ampMinImag = min(imagPert);  ampMaxImag = max(imagPert);

    lB = [ampMinReal * ones(numDipoles, 1); ampMinImag * ones(numDipoles, 1)];
    uB = [ampMaxReal * ones(numDipoles, 1); ampMaxImag * ones(numDipoles, 1)];

    
    % Create initial guess for PSO
    initialGuess = [real(dipolePerturbed.complAmpl); imag(dipolePerturbed.complAmpl)]';
    
    % Create swarm matrix by replicating the initial guess (each row is a particle)
    swarmSize = numDipoles/2;  % Swarm size is double the number of dipoles (real + imag)
    initialSwarmMatrix = repmat(initialGuess, swarmSize, 1);
    
    % Define PSO optimization settings
    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', swarmSize, ...                 % Number of particles
        'MaxIterations', 200, ...                   % Iteration limit
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
        'MaxFunctionEvaluations', 50000, ... % Maximum number of function evaluations
        'OptimalityTolerance', 1e-12, ...   % Stop if optimality conditions are met within this tolerance
        'StepTolerance', 1e-12, ...         % Stop if step size is below this threshold
        'Display', 'iter');                 % Display iteration details 
    
    optimFun = @(amp) optimization.optimFunX(amp, dip, inputData, numDipoles);
    
    [optAmps_fmincon_vec, finalError_fmincon] = fmincon(optimFun,...
        initialGuess, [], [], [], [], lB, uB, [], options_fmincon);
    
    optAmps_fmincon = optAmps_fmincon_vec(1:numDipoles).' + ...
                      1i * optAmps_fmincon_vec(numDipoles+1:end).';
    disp(['Final Error (fmincon): ', num2str(finalError_fmincon)]);
    
    
    dipoleFmincon = dipoleRef;
    dipoleFmincon.complAmpl = optAmps_fmincon;
    
%% === Far-Field Comparison: Optimized vs Reference ===
%% Step 1: Input inicialization
Ntheta = 180;
Nphi = 360;
freq = inputData.freq;
dipoleOpt = dipoleFmincon;
%% Step 1: Define Angular Grid
theta = linspace(0, pi, Ntheta).';              % Zenith angle (0 to pi)
phi   = linspace(0, 2*pi, Nphi);                % Azimuth angle (0 to 2pi)
[PHI, THETA] = meshgrid(phi, theta);            % Create 2D angular grid

%% Step 2: Convert to Cartesian Unit Vectors
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta * Nphi, 1), THETA(:), PHI(:));
r_hat = [x, y, z];                    % Flatten to N x 3 matrix

%% Step 3: Evaluate Far-Field Electric Field Vectors
fF_ref = fieldEvaluation.farFieldM2(r_hat, dipoleRef, freq);
fF_opt = fieldEvaluation.farFieldM2(r_hat, dipoleOpt, freq);

%% Step 4: Convert Cartesian to Spherical Field Components
[~, Fth_ref, Fph_ref] = utilities.transforms.cart2sphVec(x(:), y(:), z(:), fF_ref(:,1), fF_ref(:,2), fF_ref(:,3));
[~, Fth_opt, Fph_opt] = utilities.transforms.cart2sphVec(x(:), y(:), z(:), fF_opt(:,1), fF_opt(:,2), fF_opt(:,3));

%% Step 5: Reshape Components to Angular Grid
Fth_ref = reshape(Fth_ref, Ntheta, Nphi);
Fph_ref = reshape(Fph_ref, Ntheta, Nphi);
Fth_opt = reshape(Fth_opt, Ntheta, Nphi);
Fph_opt = reshape(Fph_opt, Ntheta, Nphi);

%% Step 6: Compute Total Intensity
intensity_ref = abs(Fth_ref).^2 + abs(Fph_ref).^2;
intensity_opt = abs(Fth_opt).^2 + abs(Fph_opt).^2;

%% Step 7: Normalize (Linear Scale)
intensity_ref_norm = intensity_ref / max(intensity_ref(:));
intensity_opt_norm = intensity_opt / max(intensity_opt(:));
intensity_diff = abs(intensity_opt_norm - intensity_ref_norm);

%% Step 9: Visualization
%% --- Figures: Differences ---
figure('Name','Far-Field Intensity Difference');

imagesc(phi/pi, theta/pi, intensity_diff);
xlabel('\phi/\pi'); ylabel('\theta/\pi');
title('Difference (Optimized - Reference)');
axis xy; pbaspect([2 1 1]);
colorbar;
colormap("bone");
set(gcf, 'Color', 'w');



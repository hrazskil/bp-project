%% --- 1_a Initialization ---
clc; clear; close all;

% Load dipole data from file (ensure path is correct)
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\DipoleArray.mat');

% Extract real and imaginary parts of dipole amplitudes
numDipoles = numel(dip.complAmpl);
dipoleAmpReal = real(dip.complAmpl);
dipoleAmpImag = imag(dip.complAmpl);

% Normalize the real and imaginary parts separately and store maximum values
maxAmpReal = max(abs(dipoleAmpReal));
maxAmpImag = max(abs(dipoleAmpImag));
normalizedAmpReal = dipoleAmpReal / maxAmpReal;
normalizedAmpImag = dipoleAmpImag / maxAmpImag;

% Save reference dipole amplitudes after normalization
dipoleRef = dip;
dipoleRef.complAmpl = normalizedAmpReal + 1i * normalizedAmpImag;

%% --- 1_b Perturbation of Initial Amplitudes ---
% Generate independent random perturbation factors for real and imaginary parts
realPerturbationFactor = 1 + 0.001 * randn(numDipoles, 1); % No perturbation in this case, set stddev > 0 for variation
imagPerturbationFactor = 1 + 0.001 * randn(numDipoles, 1);

% Apply perturbations separately to real and imaginary components of dipole amplitudes
perturbedAmp = real(dip.complAmpl) .* realPerturbationFactor + ...
               1i * imag(dip.complAmpl) .* imagPerturbationFactor;

% Normalize perturbed dipole amplitudes
perturbedAmpReal = real(perturbedAmp);
perturbedAmpImag = imag(perturbedAmp);
normalizedPerturbedAmpReal = perturbedAmpReal / max(abs(perturbedAmpReal));
normalizedPerturbedAmpImag = perturbedAmpImag / max(abs(perturbedAmpImag));

% Store perturbed dipole amplitudes as reference
dipolePerturbedRef.complAmpl = normalizedPerturbedAmpReal + 1i * normalizedPerturbedAmpImag;

%% --- 1_c Compute Far-Field Parameters ---
% Define physical constants using provided utility function
construct = utilities.constants.giveConstants();
omega = 2 * pi * f0List;  % Angular frequency
k = omega / construct.c0;    % Wavenumber
rFar = 1e6 / k;              % Set a large observation distance

% Get Lebedev quadrature points and weights (for integration over sphere)
Nleb = 302;                 % Number of Lebedev quadrature points (can adjust for accuracy)
[points, weights, ~] = utilities.getLebedevSphere(Nleb);
rObserved = points * rFar;   % Scale points to observation distance

% Compute the far-field pattern for the reference dipoles
fF_ref = fieldEvaluation.farField(rObserved, dipoleRef, f0List);

% Compute total radiated power for normalization
totalPower_ref = sum(sum(fF_ref .* conj(fF_ref), 2) .* weights) / (2 * construct.Z0);

%% --- 1_d Validate Objective Function Before Optimization ---
% Compute initial error using the norm of the objective function
initialError = optimization.normObjectiveFunction( ...
    dipolePerturbedRef.complAmpl, dip, f0List, points, weights, ...
    fF_ref, totalPower_ref);
disp(['Initial Test Error: ', num2str(initialError)]);

%% --- PSO Optimization for Different Numbers of Dipoles ---
numTotalDipoles = numel(dipolePerturbedRef.complAmpl);
numDipoles = round(numTotalDipoles/2);

% Randomly pick dipoles for testing
selectedIdx = randperm(numel(dip.complAmpl), numDipoles); 
dipSubset_rand.complAmpl = dipolePerturbedRef.complAmpl(selectedIdx);
dipSubset_rand.pos = dip.pos(selectedIdx, :);
dipSubset_rand.dir = dip.dir(selectedIdx, :);

% Evenly spread indices across dipoles
selectedIdx = round(linspace(1, numel(dipolePerturbedRef.complAmpl), numDipoles)); 
dipSubset_evenly.complAmpl = dipolePerturbedRef.complAmpl(selectedIdx);
dipSubset_evenly.pos = dip.pos(selectedIdx, :);
dipSubset_evenly.dir = dip.dir(selectedIdx, :);

% Sort dipoles by amplitude strength and pick the top N
[~, idxSorted] = sort(abs(dipolePerturbedRef.complAmpl), 'descend'); 
selectedIdx = idxSorted(1:numDipoles);  
dipSubset_ampl.complAmpl = dipolePerturbedRef.complAmpl(selectedIdx);
dipSubset_ampl.pos = dip.pos(selectedIdx, :);
dipSubset_ampl.dir = dip.dir(selectedIdx, :);

%% Selection Based on Error Contribution
errorContribution = zeros(numTotalDipoles, 1);

% Evaluate error contribution of each dipole
for j = 1:numTotalDipoles
    tempDip = dip;
    tempDip.complAmpl = dipolePerturbedRef.complAmpl;
    tempDip.complAmpl(j) = eps; % Set this dipole's amplitude to zero
    
    % Compute the error without this dipole
    tempError = optimization.normObjectiveFunction(tempDip.complAmpl, ...
        tempDip, f0List, points, weights, fF_ref, totalPower_ref);
    
    % Store the increase in error due to removing the dipole
    errorContribution(j) = tempError - initialError; % Compare to baseline error
end

% Sort dipoles by their impact on the error (descending order)
[~, idxSorted] = sort(errorContribution, 'descend');

% Select the top dipoles based on error contribution
selectedIdx = idxSorted(1:numDipoles);

% Subset the dipoles based on the selection
dipSubset_errorContribution.dir = dip.dir(selectedIdx, :);
dipSubset_errorContribution.pos = dip.pos(selectedIdx, :);
dipSubset_errorContribution.complAmpl = dip.complAmpl(selectedIdx);

% %% Check if no dipoles are filtered out
% disp('Checking if dipSubset is identical to dip for zero-filtering case');
% 
% % Compare the amplitudes
% maxDifference = max(abs(dipSubset.complAmpl - dipolePerturbedRef.complAmpl));
% disp(['Max difference in amplitudes: ', num2str(maxDifference)]);
% 
% % Optionally check for exact equality of the subset and original dipoles
% if maxDifference == 0
%     disp('dipSubset and dip are identical (no filtering)');
% else
%     disp('dipSubset and dip are different');
% end

%% Visualization of Dipole Selection
figure;
hold on;
grid on;
axis equal;

% Plot selected dipoles (size proportional to amplitude magnitude)
scatter3(dip.pos(selectedIdx, 1), dip.pos(selectedIdx, 2), dip.pos(selectedIdx, 3), ...
    100 * abs(dip.complAmpl(selectedIdx)) / max(abs(dip.complAmpl)), 'b', 'filled', 'DisplayName', 'Kept Dipoles');

% Plot removed dipoles (red, smaller size)
removedIdx = setdiff(1:numTotalDipoles, selectedIdx);
scatter3(dip.pos(removedIdx, 1), dip.pos(removedIdx, 2), dip.pos(removedIdx, 3), ...
    50, 'r', 'filled', 'DisplayName', 'Removed Dipoles');

xlabel('X Position');
ylabel('Y Position');
zlabel('Z Position');
title('Visualization of Dipole Selection');
legend;
view(3); % Set 3D view

%% --- Validate Objective Function after Optimization ---
% Compute the error after optimization and dipole filtering
filteredError = optimization.normObjectiveFunction( ...
    dipSubset_rand.complAmpl, dipSubset_rand, f0List, points, weights, ...
    fF_ref, totalPower_ref);
disp(['Test Error after Optimization and Dipole Filtering random: ', num2str(filteredError)]);

% Compute the error after optimization and dipole filtering
filteredError = optimization.normObjectiveFunction( ...
    dipSubset_evenly.complAmpl, dipSubset_evenly, f0List, points, weights, ...
    fF_ref, totalPower_ref);
disp(['Test Error after Optimization and Dipole Filtering evenly: ', num2str(filteredError)]);

% Compute the error after optimization and dipole filtering
filteredError = optimization.normObjectiveFunction( ...
    dipSubset_ampl.complAmpl, dipSubset_ampl, f0List, points, weights, ...
    fF_ref, totalPower_ref);
disp(['Test Error after Optimization and Dipole Filtering amplitude: ', num2str(filteredError)]);

% Compute the error after optimization and dipole filtering
filteredError = optimization.normObjectiveFunction( ...
    dipSubset_errorContribution.complAmpl, dipSubset_errorContribution, f0List, points, weights, ...
    fF_ref, totalPower_ref);
disp(['Test Error after Optimization and Dipole Filtering errorContribution: ', num2str(filteredError)]);
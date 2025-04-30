function plotFarFieldIntensityComparison(dipoleRef, dipolePso, freq, Ntheta, Nphi, refPercentiles, psoPercentiles, diffPercentiles, minRange)
% plotFarFieldIntensityComparison
% Computes and visualizes normalized far-field intensity patterns for a reference and optimized dipole configuration.
%
% INPUTS:
%   dipoleRef       - Structure containing reference dipole parameters
%   dipolePso       - Structure containing optimized dipole parameters
%   freq            - Frequency (scalar) in Hz
%   Ntheta, Nphi    - Number of samples for theta and phi angular grid
%   refPercentiles  - 2-element vector specifying color limit percentiles for reference
%   psoPercentiles  - 2-element vector specifying color limit percentiles for optimized result
%   diffPercentiles - 2-element vector for percentiles used in difference scaling
%   minRange        - Minimum range (in dB) enforced on the difference plot color scale
%
% OUTPUT:
%   Visual comparison of normalized far-field intensity in three subplots:
%     1. Reference dipole
%     2. Optimized dipole
%     3. Difference between them

%% Step 1: Define Angular Grid
theta = linspace(0, pi, Ntheta).';              % Zenith angle (0 to pi)
phi   = linspace(0, 2*pi, Nphi);                % Azimuth angle (0 to 2pi)
[PHI, THETA] = meshgrid(phi, theta);            % Create 2D angular grid

%% Step 2: Convert to Cartesian Unit Vectors
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta * Nphi, 1), THETA(:), PHI(:));
r_hat = [x, y, z];                    % Flatten to N x 3 matrix

%% Step 3: Evaluate Far-Field Electric Field Vectors
fF_ref = fieldEvaluation.farFieldM2(r_hat, dipoleRef, freq);
fF_pso = fieldEvaluation.farFieldM2(r_hat, dipolePso, freq);

%% Step 4: Convert Cartesian to Spherical Field Components
[~, Fth_ref, Fph_ref] = utilities.transforms.cart2sphVec(x(:), y(:), z(:), fF_ref(:,1), fF_ref(:,2), fF_ref(:,3));
[~, Fth_pso, Fph_pso] = utilities.transforms.cart2sphVec(x(:), y(:), z(:), fF_pso(:,1), fF_pso(:,2), fF_pso(:,3));

%% Step 5: Reshape Components to Angular Grid
Fth_ref = reshape(Fth_ref, Ntheta, Nphi);
Fph_ref = reshape(Fph_ref, Ntheta, Nphi);
Fth_pso = reshape(Fth_pso, Ntheta, Nphi);
Fph_pso = reshape(Fph_pso, Ntheta, Nphi);

%% Step 6: Compute Total Intensity
intensity_ref = abs(Fth_ref).^2 + abs(Fph_ref).^2;
intensity_pso = abs(Fth_pso).^2 + abs(Fph_pso).^2;

%% Step 7: Convert to dB Scale and Normalize
intensity_ref_dB = 10 * log10(intensity_ref + eps);
intensity_pso_dB = 10 * log10(intensity_pso + eps);

intensity_ref_dB_norm = intensity_ref_dB - max(intensity_ref_dB(:));
intensity_pso_dB_norm = intensity_pso_dB - max(intensity_pso_dB(:));
intensity_diff_dB = intensity_pso_dB_norm - intensity_ref_dB_norm;

%% Step 8: Compute Smart Color Axis Limits
clim_ref = [-max(abs(prctile(intensity_ref_dB_norm(:), refPercentiles))) 0];
clim_pso = [-max(abs(prctile(intensity_pso_dB_norm(:), psoPercentiles))) 0];

vmin = prctile(intensity_diff_dB(:), diffPercentiles(1));
vmax = prctile(intensity_diff_dB(:), diffPercentiles(2));
vmax_abs = max([abs(vmin), abs(vmax), minRange]);
clim_diff = [-vmax_abs, vmax_abs];

%% Step 9: Visualization
figure('Name','Far-Field Intensity Comparison (Normalized)', 'Position', [100, 100, 1400, 500]);

subplot(1,3,1);
imagesc(rad2deg(phi), rad2deg(theta), intensity_ref_dB_norm);
xlabel('\phi [\circ]'); ylabel('\theta [\circ]');
title('Reference Intensity (Normalized dB)');
colorbar; axis xy; clim(clim_ref);

subplot(1,3,2);
imagesc(rad2deg(phi), rad2deg(theta), intensity_pso_dB_norm);
xlabel('\phi [\circ]'); ylabel('\theta [\circ]');
title('Optimized Intensity (Normalized dB)');
colorbar; axis xy; clim(clim_pso);

subplot(1,3,3);
imagesc(rad2deg(phi), rad2deg(theta), intensity_diff_dB);
xlabel('\phi [\circ]'); ylabel('\theta [\circ]');
title('Difference (Optimized - Reference) [dB]');
colorbar; axis xy; clim(clim_diff);
end


function plotFarFieldIntensityComparison(dipoleRef, dipoleOpt, freq, Ntheta, Nphi)
% plotFarFieldIntensityComparison
% Computes and visualizes normalized far-field intensity patterns for a reference and optimized dipole configuration.
%
% INPUTS:
%   dipoleRef       - Structure containing reference dipole parameters
%   dipoleOpt       - Structure containing optimized dipole parameters
%   freq            - Frequency (scalar) in Hz
%   Ntheta, Nphi    - Number of samples for theta and phi angular grid
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
% figure('Name','Far-Field Intensity Comparison (Normalized)');
% colormap("bone");  % Zvol intenzitní mapu pro celý obrázek
% 
% % --- Plot 1: Reference ---
% subplot(3,1,1);
% imagesc(phi/pi, theta/pi, intensity_ref_norm);
% xlabel('\phi/\pi'); ylabel('\theta/\pi');
% title('Reference Intensity');
% colorbar; axis xy;
% % Capture the color axis limits from the first plot
% refCaxis = clim;
% pbaspect([2 1 1]);       % Set aspect ratio
% 
% % --- Plot 2: Optimized ---
% subplot(3,1,2);
% imagesc(phi/pi, theta/pi, intensity_opt_norm);
% xlabel('\phi/\pi'); ylabel('\theta/\pi');
% title('Optimized Intensity');
% colorbar; axis xy;
% clim(refCaxis);  % Apply reference color scale
% pbaspect([2 1 1]);       % Set aspect ratio

% %% --- Figure 1: Reference and Optimized (Shared Colorbar & CLim) ---
% figure('Name','Far-Field Intensity: Reference vs Optimized');
% 
% tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact'); % Only 2 tiles
% 
% % --- Reference ---
% ax1 = nexttile;
% imagesc(phi/pi, theta/pi, intensity_ref_norm);
% xlabel('\phi/\pi'); ylabel('\theta/\pi');
% title('Reference Intensity');
% axis xy; pbaspect([2 1 1]);
% refCaxis = clim;  % Save caxis
% 
% % --- Optimized ---
% ax2 = nexttile;
% imagesc(phi/pi, theta/pi, intensity_opt_norm);
% xlabel('\phi/\pi'); ylabel('\theta/\pi');
% title('Optimized Intensity');
% axis xy; pbaspect([2 1 1]);
% clim(refCaxis);  % Match color scale
% 
% % --- Shared Colorbar ---
% cb = colorbar(ax2, 'Location', 'eastoutside');
% cb.Label.String = 'Normalized Intensity';
% 
% % Match clim
% clim(ax1, refCaxis);
% clim(ax2, refCaxis);
% 
% colormap("bone");
% set(gcf, 'Color', 'w');

%% --- Figure 2: Difference ---
figure('Name','Far-Field Intensity Difference');

imagesc(phi/pi, theta/pi, intensity_diff);
xlabel('\phi/\pi'); ylabel('\theta/\pi');
title('Difference (Optimized - Reference)');
axis xy; pbaspect([2 1 1]);
colorbar;
colormap("bone");
set(gcf, 'Color', 'w');
end


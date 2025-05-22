 %% --- Initialization ---
    clc; clear; close all;

load('results/test_structure_2_N_dip_scaling.mat', 'results');
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat');

Ntheta = 90;
Nphi = 180;
phi = linspace(0, 2*pi, Nphi);
theta = linspace(0, pi, Ntheta).';

NdipList = [results.Ndip];

intensity_diff_all = cell(1, numel(results));
max_vals = zeros(1, numel(results));

frequency = f0List;

% === Precompute reference far-field once ===
[PHI, THETA] = meshgrid(phi, theta);
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta*Nphi,1), THETA(:), PHI(:));
r_hat = [x, y, z];

% Evaluate reference far-field
fF_ref = fieldEvaluation.farFieldM2(r_hat, dip, frequency);
[~, Fth_ref, Fph_ref] = utilities.transforms.cart2sphVec(x, y, z, ...
    fF_ref(:,1), fF_ref(:,2), fF_ref(:,3));

% Compute and reshape reference intensity
I_ref = abs(Fth_ref).^2 + abs(Fph_ref).^2;
I_ref = reshape(I_ref, Ntheta, Nphi);
I_ref_norm = I_ref / max(I_ref(:));  % normalize once

% === Loop over results ===
for i = 1:numel(results)
    % Rebuild optimized dipole
    dipOpt.complAmpl = results(i).optAmps_fmincon;
    dipOpt.pos = dip.pos(results(i).selectedIdx, :);
    dipOpt.dir = dip.dir(results(i).selectedIdx, :);

    % Evaluate optimized far-field
    fF_opt = fieldEvaluation.farFieldM2(r_hat, dipOpt, frequency);
    [~, Fth_opt, Fph_opt] = utilities.transforms.cart2sphVec(x, y, z, ...
        fF_opt(:,1), fF_opt(:,2), fF_opt(:,3));

    I_opt = abs(Fth_opt).^2 + abs(Fph_opt).^2;
    I_opt = reshape(I_opt, Ntheta, Nphi);
    I_opt_norm = I_opt / max(I_opt(:));

    % Difference in normalized intensity
    diff = abs(I_opt_norm - I_ref_norm);
    intensity_diff_all{i} = diff;
    max_vals(i) = max(diff(:));
end

% Use max of last case or global max for consistency
clim = [0, max(max_vals)];

% Plotting
nCols = 2;
nRows = ceil(numel(results)/nCols);
figure('Name', 'Far-Field Difference Comparison', 'Color', 'w');

for i = 1:numel(results)
    subplot(nRows, nCols, i);
    imagesc(phi/pi, theta/pi, intensity_diff_all{i});
    title(sprintf('%d Dipoles', results(i).Ndip));
    xlabel('\phi/\pi'); ylabel('\theta/\pi');
    axis xy; pbaspect([2 1 1]);
    caxis(clim);
end

% Shared colorbar
cb = colorbar('Position', [0.93 0.1 0.02 0.8]);
cb.Label.String = 'Normalized Difference';
colormap('bone');

 %% --- Initialization ---
    clc; clear; close all;

load('results/test_structure_2_N_dip_scaling.mat', 'results');
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat');

Ntheta = 180;
Nphi = 360;
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

% === LaTeX-compatible defaults ===
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% === PLOTTING ===
% Layout constants
nCols = 4;
subplotWidth = 0.14;
subplotHeight = 0.8;
hSpacing = 0.015;
leftMargin = 0.02;
bottomMargin = 0.15;
cbWidth = 0.015;
cbHeight = 0.75;

% Overall figure
figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 6]);

% Compute global clim for consistency
clim = [0, max(cellfun(@(A) max(A(:)), intensity_diff_all(1:4)))];

% Plot each subplot
for j = 1:nCols
    left = leftMargin + (j-1)*(subplotWidth + hSpacing);
    ax = axes('Position', [left, bottomMargin, subplotWidth, subplotHeight]);

    imagesc(theta/pi, phi/pi, intensity_diff_all{j});
    axis xy;
    pbaspect([1 2 1]);
    colormap('bone');
    caxis(clim);

    if j == 1
        xlabel('$\theta/\pi$', 'Interpreter', 'latex');
        ylabel('$\phi/\pi$', 'Interpreter', 'latex');
    else
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
    end
end

% Set axes units for TikZ compatibility
set(findall(gcf, 'Type', 'axes'), 'Units', 'centimeters');

% Add colorbar
cbLeft = left + subplotWidth + 0.01;
cbBottom = bottomMargin + 0.02;
cb = colorbar('Position', [cbLeft, cbBottom, cbWidth, cbHeight]);
cb.Label.String = 'Normalized Difference';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;

% === EXPORT ===
matlab2tikz('FF_I_Ref_vs_Opt_Diff_dip_1.tikz', ...
    'width', '\figurewidth', ...
    'height', '\figureheight', ...
    'showInfo', false, ...
    'extraAxisOptions', {
        'tick label style={font=\scriptsize}', ...
        'label style={font=\normalsize}', ...
        'title style={font=\normalsize}', ...
        'colormap name=bone' ...
    });

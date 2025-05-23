 %% --- Initialization for test_structure_2 ---
    clc; clear; close all;

load('test_structure_2_N_obs_scaling.mat', 'results');
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat');


%% Grid settings
Ntheta = 180;
Nphi = 360;
theta = linspace(0, pi, Ntheta).';
phi   = linspace(0, 2*pi, Nphi);
[PHI, THETA] = meshgrid(phi, theta);
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta * Nphi, 1), THETA(:), PHI(:));
r_hat = [x, y, z];

%% Reference dipole (use same as optimization)
dipoleFull = dip;  % from your DipoleArray.mat
freq = f0List;  % assumes freq is the same
fF_ref = fieldEvaluation.farFieldM2(r_hat, dipoleFull, freq);
[~, Fth_ref, Fph_ref] = utilities.transforms.cart2sphVec(x, y, z, fF_ref(:,1), fF_ref(:,2), fF_ref(:,3));
I_ref = abs(Fth_ref).^2 + abs(Fph_ref).^2;
I_ref = reshape(I_ref, Ntheta, Nphi);
I_ref_norm = I_ref / max(I_ref(:));

%% Evaluate each result
range_vals = zeros(1, numel(results));
intensity_diff_all = cell(1, numel(results));

for i = 1:numel(results)
    dipOpt.complAmpl = results(i).optAmps_fmincon;
    dipOpt.pos = dip.pos;
    dipOpt.dir = dip.dir;

    fF_opt = fieldEvaluation.farFieldM2(r_hat, dipOpt, freq);
    [~, Fth_opt, Fph_opt] = utilities.transforms.cart2sphVec(x, y, z, ...
        fF_opt(:,1), fF_opt(:,2), fF_opt(:,3));

    I_opt = abs(Fth_opt).^2 + abs(Fph_opt).^2;
    I_opt = reshape(I_opt, Ntheta, Nphi);
    I_opt_norm = I_opt / max(I_opt(:));

    diff = abs(I_opt_norm - I_ref_norm);
    intensity_diff_all{i} = diff;
    range_vals(i) = max(diff(:)) - min(diff(:));
end
%%
% Set LaTeX-compatible interpreters globally
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Set parameters
nCols = numel(results);
clim = [0, max(range_vals)];

% Layout parameters
subplotWidth = 0.14;
subplotHeight = 0.8;
hSpacing = 0.015;
leftMargin = 0.02;
bottomMargin = 0.15;

% Create figure
figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 6]);

for i = 1:nCols
    % Compute position of each subplot
    left = leftMargin + (i-1)*(subplotWidth + hSpacing);
    ax = axes('Position', [left, bottomMargin, subplotWidth, subplotHeight]);

    % Plot data
    imagesc(theta/pi, phi/pi, intensity_diff_all{i});
    axis xy;
    pbaspect([1 2 1]);
    colormap('bone');
    caxis(clim);

    % Labels only on first
    if i == 1
        xlabel('$\theta/\pi$', 'Interpreter', 'latex');
        ylabel('$\phi/\pi$', 'Interpreter', 'latex');
    else
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
    end
end

% Ensure all axes are in physical units (not 'normalized')
set(findall(gcf,'Type','axes'), 'Units', 'centimeters');

% Add shared colorbar manually
cbWidth = 0.015;
cbHeight = 0.75;
cbLeft = left + subplotWidth + 0.01;  % after last plot
cbBottom = bottomMargin + 0.025;

cb = colorbar('Position', [cbLeft, cbBottom, cbWidth, cbHeight]);
cb.Label.String = 'Normalized Difference';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;

% --- Export to TikZ ---
matlab2tikz('FF_I_Ref_vs_Per_Diff_obs.tikz', ...
    'width', '\figurewidth', ...
    'height', '\figureheight', ...
    'showInfo', false, ...
    'extraAxisOptions', {
        'tick label style={font=\scriptsize}', ...
        'label style={font=\normalsize}', ...
        'title style={font=\normalsize}', ...
        'colormap name=bone' ...
    });
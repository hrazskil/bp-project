 %% --- Initialization for test_structure_2 ---
    clc; clear; close all;

load('test_structure_2_N_obs_scaling.mat', 'results');
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat');


%% Grid settings
Ntheta = 90;
Nphi = 180;
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

%% Plotting
clim = [0, max(range_vals)];
NobsList = [results.Nobs];
nCols = 3;
nRows = ceil(numel(results)/nCols);

figure('Name', 'Far-Field Comparison: Varying Nobs', 'Color', 'w');
for i = 1:numel(results)
    subplot(nRows, nCols, i);
    imagesc(phi/pi, theta/pi, intensity_diff_all{i});
    title(sprintf('%d Obs Points', NobsList(i)));
    xlabel('\phi/\pi'); ylabel('\theta/\pi');
    axis xy; pbaspect([2 1 1]);
    caxis(clim);
end

% Shared colorbar
cb = colorbar('Position', [0.93 0.1 0.02 0.8]);
cb.Label.String = 'Normalized Difference';
colormap('bone');

% Split into groups of 4
groupSize = 4;
numGroups = ceil(numel(results) / groupSize);

for g = 1:numGroups
    figure('Name', sprintf('Far-Field Difference Comparison: Group %d', g), 'Color', 'w');

    % Determine indices for this group
    idxStart = (g-1)*groupSize + 1;
    idxEnd = min(g*groupSize, numel(results));
    groupIdx = idxStart:idxEnd;
    
    for j = 1:length(groupIdx)
        subplot(2, 2, j);  % 2 rows x 2 columns for 4 plots
        i = groupIdx(j);

        imagesc(phi/pi, theta/pi, intensity_diff_all{i});
        title(sprintf('%d Dipoles', results(i).Ndip));
        xlabel('\phi/\pi'); ylabel('\theta/\pi');
        axis xy; pbaspect([2 1 1]);
        caxis(clim);
    end

    % Shared colorbar to the right of subplot (approximate position)
    cb = colorbar('Position', [0.91 0.1 0.02 0.8]);
    cb.Label.String = 'Normalized Difference';
    colormap('bone');
end
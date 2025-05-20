function plotFarFieldComparison(dipoleRef, dipoleSim, freq, Ntheta, Nphi)
% plotFarFieldComparison
% Compares far-field radiation patterns.
%
% INPUTS:
%   dipoleRef    - Reference dipole structure
%   dipolePso    - Optimized dipole structure
%   freq         - Frequency in Hz
%   Ntheta, Nphi - Number of angular samples (theta, phi)
%
% OUTPUT:
%   Contour plot visualizig relative error.

%% Step 1: Generate Angular Grid
theta = linspace(0, pi, Ntheta);
phi   = linspace(0, 2*pi, Nphi);
[PHI, THETA] = meshgrid(phi, theta);

%% Step 2: Convert to Cartesian Unit Vectors
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta * Nphi, 1), THETA(:), PHI(:));

%% Step 3: Evaluate Far-Field Vectors
fF_ref = fieldEvaluation.farFieldM2([x, y, z], dipoleRef, freq);
fF_sim = fieldEvaluation.farFieldM2([x, y, z], dipoleSim, freq);

%% Step 5: Error computation
% Error = abs(fF_sim - fF_ref)/(max(abs(fF_sim)));

Error = utilities.rowNorm(abs(fF_sim - fF_ref)) / max(utilities.rowNorm(fF_ref));

%% Step 6: Plot Comparison of Each Component
ErrorMap = reshape(Error, Ntheta, Nphi);

    figure;
    contourf(PHI/pi, THETA/pi, ErrorMap, 50, 'LineColor', 'none');
    xlabel('\phi/\pi'); ylabel('\theta/\pi');
    title(['Error between simulated and reference Farfield']);
    colorbar;
end
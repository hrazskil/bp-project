function plotFarFieldComponentComparison(dipoleRef, dipolePso, freq, Ntheta, Nphi)
% plotFarFieldComponentComparison
% Compares vector components (F_theta and F_phi) of far-field radiation patterns.
%
% INPUTS:
%   dipoleRef    - Reference dipole structure
%   dipolePso    - Optimized dipole structure
%   freq         - Frequency in Hz
%   Ntheta, Nphi - Number of angular samples (theta, phi)
%
% OUTPUT:
%   Series of contour plots visualizing real and imaginary parts of F_theta and F_phi
%   for both dipoles on a normalized angular grid.

%% Step 1: Generate Angular Grid
theta = linspace(0, pi, Ntheta);
phi   = linspace(0, 2*pi, Nphi);
[PHI, THETA] = meshgrid(phi, theta);

%% Step 2: Convert to Cartesian Unit Vectors
[x, y, z] = utilities.transforms.sph2cartCoor(ones(Ntheta * Nphi, 1), THETA(:), PHI(:));
r_hat = [x, y, z];

%% Step 3: Evaluate Far-Field Vectors
fF_ref = fieldEvaluation.farFieldM2(r_hat, dipoleRef, freq);
fF_pso = fieldEvaluation.farFieldM2(r_hat, dipolePso, freq);

%% Step 4: Transform to Spherical Field Components
[~, Fth_ref, Fph_ref] = utilities.transforms.cart2sphVec(x, y, z, fF_ref(:,1), fF_ref(:,2), fF_ref(:,3));
[~, Fth_pso, Fph_pso] = utilities.transforms.cart2sphVec(x, y, z, fF_pso(:,1), fF_pso(:,2), fF_pso(:,3));

%% Step 5: Reshape to Angular Grid
Fth_ref = reshape(Fth_ref, Ntheta, Nphi);
Fph_ref = reshape(Fph_ref, Ntheta, Nphi);
Fth_pso = reshape(Fth_pso, Ntheta, Nphi);
Fph_pso = reshape(Fph_pso, Ntheta, Nphi);

%% Step 6: Plot Comparison of Each Component
components = {'F_\theta', 'F_\phi'};
fields_ref = {Fth_ref, Fph_ref};
fields_pso = {Fth_pso, Fph_pso};

for i = 1:2
    figure;
    subplot(1,2,1);
    contourf(PHI/pi, THETA/pi, real(fields_ref{i}));
    xlabel('\theta/\pi'); ylabel('\phi/\pi');
    title(['Reference Re[', components{i}, ']']);
    colorbar;

    subplot(1,2,2);
    contourf(PHI/pi, THETA/pi, real(fields_pso{i}));
    xlabel('\theta/\pi'); ylabel('\phi/\pi');
    title(['Optimized Re[', components{i}, ']']);
    colorbar;

    figure;
    subplot(1,2,1);
    contourf(PHI/pi, THETA/pi, imag(fields_ref{i}));
    xlabel('\theta/\pi'); ylabel('\phi/\pi');
    title(['Reference Im[', components{i}, ']']);
    colorbar;

    subplot(1,2,2);
    contourf(PHI/pi, THETA/pi, imag(fields_pso{i}));
    xlabel('\theta/\pi'); ylabel('\phi/\pi');
    title(['Optimized Im[', components{i}, ']']);
    colorbar;
end
end

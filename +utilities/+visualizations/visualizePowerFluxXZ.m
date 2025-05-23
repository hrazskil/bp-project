function visualizePowerFluxXZ(dip, f, xLim, zLim, nPoints)
% Visualize vector power flux (Poynting vector) in xz-plane (y=0)

% Step 1: Define grid
x = linspace(xLim(1), xLim(2), nPoints);
z = linspace(zLim(1), zLim(2), nPoints);
[X, Z] = meshgrid(x, z);
Y = zeros(size(X)); % y = 0 plane
rObserved = [X(:), Y(:), Z(:)];

% Step 2: Compute E and H fields
eF = fieldEvaluation.eleFieldM2(rObserved, dip, f);
mF = fieldEvaluation.magFieldM2(rObserved, dip, f);

% Step 3: Compute Poynting vector and extract S_x
S = fieldEvaluation.powerPoynting(eF, mF); % (nObs × 3)
Sx = reshape(S(:,1), size(X));

% Step 4: Plot with imagesc
figure;
imagesc(x, z, Sx);
axis xy;
xlabel('x [m]');
ylabel('z [m]');
title('S_x Mean (W/m^2) — Power Flux in x-Direction');
colorbar;
colormap('hot');  % or 'parula', 'bone', etc.
end
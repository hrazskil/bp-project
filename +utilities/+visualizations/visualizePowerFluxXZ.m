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

% Step 3: Compute Poynting vector
S = fieldEvaluation.powerPoynting(eF, mF); % (nObs × 3)

% Step 4: Reshape for plotting
Sx = reshape(S(:,1), size(X));
Sz = reshape(S(:,3), size(Z));

% Step 5: Plot vector field
figure;
quiver(X, Z, Sx, Sz, 'k');
axis equal;
xlabel('x [m]');
ylabel('z [m]');
title('Power Flux (Poynting Vector) in xz-Plane');
end
function PowerFluxXZ(dip, f, xLim, zLim, nPoints)
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

% Step 4: Set cutoff values to NaN
Sx(Sx < -250 | Sx > 250) = NaN;

% Step 5: Plot with contourf
levels = linspace(-250, 250, 50);  % Define the contour levels

%Filled contour plot (no edge lines)
figure;

% Step 1: Filled contours (without edge lines)
contourf(x, z, Sx, levels, 'LineWidth', 1);
colormap(parula);       % or customMap
colorbar;
axis xy;
xlabel('x [m]');
ylabel('z [m]');
title('S_x Mean (W/m^2)');
end
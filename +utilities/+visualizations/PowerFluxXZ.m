function PowerFluxXZ(dip, f, xLim, zLim, xPoints, zPoints)
% Visualize vector power flux (Poynting vector) in xz-plane (y=0)

% Step 1: Define grid
x = linspace(xLim(1), xLim(2), xPoints);
z = linspace(zLim(1), zLim(2), zPoints);
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
    
%Filled contour plot (no edge lines)
figure;

% Step 1: Filled contours (without edge lines)
contourf(x,z,Sx,[linspace(-40,-1,30),linspace(10,250,30)]);
colormap(parula);       % or customMap
colorbar;
clim([-250 250]);
% axis xy;
xlabel('x [m]');
ylabel('z [m]');
title('S_x Mean (W/m^2)');

% Remove everything: axis, ticks, labels, whitespace
axis off;
axis tight;
set(gca, 'Position', [0 0 1 1]);       % Fill entire figure
set(gca, 'LooseInset', [0 0 0 0]);     % No margin

% Export as clean PNG
exportgraphics(gcf, 'Sx_ref.png', 'Resolution', 300);

end
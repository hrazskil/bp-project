clc;
clear;

%% Initialization for alpha testing transformation functions
[x, y, z] = ndgrid([-1, 0, 1]);
rObserved = [x(:), y(:), z(:)];

%% inicialization for beta testing transformation functions
% nObs=100;
% rObserved=rand(nObs,3);
% 
% rObsZero=rObserved.*(1./(rObserved(:,1)+rObserved(:,2)+rObserved(:,3)));

%% inicialization for gamma testing transformation functions
% nObs=100;
% Define the number of theta and phi points
nTh = 10; 
nPh = 49;

% Create linearly spaced vectors for theta and phi
theta = linspace(0, pi, nTh); 
phi = linspace(0, 2*pi, nPh);

% Create a meshgrid for theta and phi
[Theta, Phi] = meshgrid(theta, phi);

% Reshape Theta and Phi into column vectors
Theta = reshape(Theta, [nTh*nPh, 1]);
Phi = reshape(Phi, [nTh*nPh, 1]);

% Convert spherical coordinates to Cartesian coordinates
[x, y, z] = utilities.transforms.sph2cartCoor(ones(nTh*nPh, 1), Theta, Phi);
rObserved = [x, y, z]; % Combine Cartesian coordinates into a matrix
%% Testing transformation functions
[sphRObs(:,1), sphRObs(:,2), sphRObs(:,3)] = ...
    utilities.transforms.cart2sphCoor(rObserved(:,1), rObserved(:,2), rObserved(:,3));

F = ones(size(rObserved,1),3);


[r0, theta0, phi0] = ...
    utilities.transforms.cart2sphVec(rObserved(:,1), rObserved(:,2), rObserved(:,3),F(:,1),F(:,2),F(:,3));

[x0, y0, z0] = ...
    utilities.transforms.sph2cartVec(sphRObs(:,2), sphRObs(:,3),r0, theta0, phi0);

%% Visualization of Original Coordinates
sRObs = size(rObserved, 1);

figure;
quiver3(zeros(sRObs, 1), zeros(sRObs, 1), zeros(sRObs, 1), ...
    rObserved(:,1), rObserved(:,2), rObserved(:,3), 'filled', 'k');
hold on;

% Unit vectors
quiver3(rObserved(:,1), rObserved(:,2), rObserved(:,3), ...
    theta0(:,1), theta0(:,2), theta0(:,3), 'b', 'LineWidth', 1);

% Labels and legend
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Scatter Plot of Observed and Transformed Coordinates');
legend('Observed', 'Unit Vectors');
grid on;
hold off;

%% Visualization of Transformed Coordinates

figure;
quiver3(zeros(sRObs, 1), zeros(sRObs, 1), zeros(sRObs, 1), ...
    rObserved(:,1), rObserved(:,2), rObserved(:,3), 'filled', 'k');
hold on;

% Unit vectors
quiver3(rObserved(:,1), rObserved(:,2), rObserved(:,3), ...
    theta02(:,1), theta02(:,2), theta02(:,3), 'b', 'LineWidth', 1);

% Labels and legend
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Scatter Plot of Observed and Transformed Coordinates');
legend('Observed', 'Unit Vectors');
grid on;
hold off;
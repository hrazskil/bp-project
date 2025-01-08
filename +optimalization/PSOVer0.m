% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%% Getting the field ready
% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\test\halfwaveDipole.mat');

% Define the number of theta and phi points
nTh = 50; 
nPh = 99;

% Create linearly spaced vectors for theta and phi
theta = linspace(0, pi, nTh); 
phi = linspace(0, 2*pi, nPh);

% Create a meshgrid for theta and phi
[Theta, Phi] = meshgrid(theta, phi);

% Reshape Theta and Phi into column vectors
Theta = reshape(Theta, [nTh*nPh, 1]);
Phi = reshape(Phi, [nTh*nPh, 1]);

% Convert spherical coordinates to Cartesian coordinates
[x, y, z] = utilities.transforms.sph2cartCoor(ones(nTh*nPh, 1)*, Theta, Phi);
rObserved = [x, y, z]; % Combine Cartesian coordinates into a matrix

% Evaluate the far field based on the observed coordinates
[fF] = fieldEvaluation.farField(rObserved, dip, f0List);


%% --- Initialization for BTS data without DipoleArray.mat ---
clc; clear; close all;

%% Load saved far-field data
load('tests/test_structure_2/dataVertical.tsv');   % vertical far-field data
load('tests/test_structure_2/dataHorizontal.tsv'); % horizontal far-field data

%% Placeholder: manually define dipole configuration and frequency
% NOTE: You need to define 'dip' structure elsewhere before optimization.
% This includes:
%   dip.pos        (N×3 matrix of positions),
%   dip.dir        (N×3 matrix of unit orientation vectors),
%   dip.complAmpl  (N×1 complex vector of amplitudes)
f0List = 2e9;  % Example frequency in Hz (adjust as appropriate)

%% Number of observation points
Nbeta = size(dataHorizontal(:,2),1);  % number of azimuthal observations
Nphi = size(dataVertical(:,2),1);     % number of elevation observations

% Reconstruct Cartesian observation points from angles
x1 = cos(dataVertical(:,1));
y1 = 0*x1;
z1 = sin(dataVertical(:,1));

x2 = cos(dataHorizontal(:,1));
y2 = sin(dataHorizontal(:,1));
z2 = 0*x2;

vertical_points = [x1, y1, z1];
horizontal_points = [x2, y2, z2];

% Create input data structure for optimization
inputData.horizontal.rad = dataHorizontal(:,2);
inputData.horizontal.points = horizontal_points;
inputData.horizontal.weights = ones(Nbeta, 1);
inputData.horizontal.totalPower = sum(inputData.horizontal.rad);

inputData.vertical.rad = dataVertical(:,2);
inputData.vertical.points = vertical_points;
inputData.vertical.weights = ones(Nphi, 1);
inputData.vertical.totalPower = sum(inputData.vertical.rad);

inputData.freq = f0List;

% generating dipole array



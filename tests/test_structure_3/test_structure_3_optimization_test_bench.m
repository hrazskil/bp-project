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

%% generating dipole array


% Parameters
numPairs = 4;
dipoleLength = 0.05;                 % Total dipole length [m]
numDipoles = 50;
offsetStep = 0.05;                   % Offset between dipole pairs [m]
rotationY = pi/4;                    % Rotation for y-direction dipoles
rotationX = 3*pi/4;                  % Rotation for x-direction dipoles


% Optional: define excitation amplitudes
amplitudes = ones(numDipoles, 1);  % Uniform excitation

% Initialize empty structure for full array
dip.pos = [];
dip.dir = [];
dip.complAmpl = [];

% Generate crossed dipole pairs
for i = 1:numPairs
    % Vertical z-offset for stacking
    zOffset = (i - 1) * offsetStep;
    posOffset = [0, 0, zOffset];

    % Generate dipole rotated for Y orientation
    dipY = geometry.halfwaveDipoleArray(dipoleLength, numDipoles, rotationY, posOffset, amplitudes);

    % Generate dipole rotated for X orientation
    dipX = geometry.halfwaveDipoleArray(dipoleLength, numDipoles, rotationX, posOffset, amplitudes);

    % Append to full structure
    dip.pos        = [dip.pos; dipY.pos; dipX.pos];
    dip.dir        = [dip.dir; dipY.dir; dipX.dir];
    dip.complAmpl  = [dip.complAmpl; dipY.complAmpl; dipX.complAmpl];
end


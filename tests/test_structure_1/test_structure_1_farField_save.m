%% --- 1_a Initialization ---
clc; clear; close all;

%% Load dipole configuration
% Loads dipole parameters including positions and amplitudes
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\halfwaveDipole.mat')

%% --- 1_c Compute Far-Field Parameters (Uniform Grid) ---
% Define number of sampling points for two angular dimensions
Nbeta = 360;  % Number of points for vertical plane (elevation)
Nphi = 360;   % Number of points for horizontal plane (azimuth)

% Angular values for elevation and azimuth
beta_vals = linspace(0, 2*pi, Nbeta).';
phi_vals = linspace(0, 2*pi, Nphi).';

% Generate unit vectors for vertical observation plane
x1 = cos(beta_vals);
y1 = 0*x1;
z1 = sin(beta_vals);

% Generate unit vectors for horizontal observation plane
x2 = cos(phi_vals);
y2 = sin(phi_vals);
z2 = 0*x2;

% Combine coordinates into point matrices
vertical_points = [x1(:,1), y1(:,1), z1(:,1)];
horizontal_points = [x2(:,1), y2(:,1), z2(:,1)];

% Evaluate far-field radiation pattern for both observation planes
farField_plane1 = fieldEvaluation.farFieldM2(vertical_points, dip, f0List);  
farField_plane2 = fieldEvaluation.farFieldM2(horizontal_points, dip, f0List);  

% Compute radiated power (squared magnitude of far-field)
rad_Vertical = sum(abs(farField_plane1).^2, 2);
rad_Horizontal = sum(abs(farField_plane2).^2, 2);

% Save computed data for further processing
data = [beta_vals, rad_Vertical];
fid = fopen('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\dataVertical.tsv', 'w');
fprintf(fid, '%.16e\t%.16e\n', data.');
fclose(fid);

data = [phi_vals, rad_Horizontal];
fid = fopen('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_1\dataHorizontal.tsv', 'w');
fprintf(fid, '%.16e\t%.16e\n', data.');
fclose(fid);
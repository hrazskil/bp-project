%% --- test_structure_2: Save Far-Field TSV Files ---
clc; clear; close all;

%% Load dipole array structure and frequency
load('tests/test_structure_2/DipoleArray.mat');  % loads 'dip' and 'f0List'
f0 = f0List;

%% Generate far-field for vertical and horizontal planes
N_points = 90;
angles = linspace(0, 2*pi*(N_points-1)/N_points, N_points).';

horizontal_points = [cos(angles), sin(angles), zeros(N_points,1)];
vertical_points   = [cos(angles), zeros(N_points,1), sin(angles)];

E_horiz = fieldEvaluation.farFieldM2(horizontal_points, dip, f0);
I_horiz = fieldEvaluation.powerDensityFar(E_horiz);

E_vert  = fieldEvaluation.farFieldM2(vertical_points, dip, f0);
I_vert  = fieldEvaluation.powerDensityFar(E_vert);

%% Save to .tsv files
% Horizontal
dataH = [angles, I_horiz];
fidH = fopen('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\dataHorizontal.tsv', 'w');
fprintf(fidH, '%.16e\t%.16e\n', dataH.');
fclose(fidH);

% Vertical
dataV = [angles, I_vert];
fidV = fopen('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\dataVertical.tsv', 'w');
fprintf(fidV, '%.16e\t%.16e\n', dataV.');
fclose(fidV);

disp('Far-field .tsv files saved successfully for test_structure_2.');

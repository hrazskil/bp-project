%% --- test_structure_2: Save Far-Field TSV Files with Lebedev Points ---
clc; clear; close all;

%% Load dipole array structure and frequency
load('tests/test_structure_2/DipoleArray.mat');  % loads 'dip' and 'f0List'
f0 = f0List;

%% Define the Lebedev degrees (number of observation points)
lebedevDegrees = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810];

%% Generate far-field data for each Lebedev degree
for degIdx = 1:length(lebedevDegrees)
    N_points = lebedevDegrees(degIdx)/2;  % Set the number of points based on the current Lebedev degree/2 so as to maintain same number of points for al planes against lebedev
    angles = linspace(0, 2*pi*(N_points-1)/N_points, N_points).';  % Evenly spaced angles for the cut

    % Define the horizontal (xy-plane) and vertical (xz-plane) cuts
    horizontal_points = [cos(angles), sin(angles), zeros(N_points,1)];
    vertical_points   = [cos(angles), zeros(N_points,1), sin(angles)];

    % Calculate the far-field for horizontal and vertical cuts
    E_horiz = fieldEvaluation.farFieldM2(horizontal_points, dip, f0);
    I_horiz = fieldEvaluation.powerDensityFar(E_horiz);

    E_vert  = fieldEvaluation.farFieldM2(vertical_points, dip, f0);
    I_vert  = fieldEvaluation.powerDensityFar(E_vert);

    %% Save the data to .tsv files with the new path
    % Horizontal (xy-plane)
    dataH = [angles, I_horiz];
    filenameH = sprintf('C:\\Users\\kilia\\Plocha\\gitHub\\bp-project\\tests\\test_structure_2\\dataHV_N_obs_points\\dataHV_%d_obs_points_Horizontal.tsv', N_points);
    fidH = fopen(filenameH, 'w');
    fprintf(fidH, '%.16e\t%.16e\n', dataH.');
    fclose(fidH);

    % Vertical (xz-plane)
    dataV = [angles, I_vert];
    filenameV = sprintf('C:\\Users\\kilia\\Plocha\\gitHub\\bp-project\\tests\\test_structure_2\\dataHV_N_obs_points\\dataHV_%d_obs_points_Vertical.tsv', N_points);
    fidV = fopen(filenameV, 'w');
    fprintf(fidV, '%.16e\t%.16e\n', dataV.');
    fclose(fidV);

    % Display status message
    disp(['Saved far-field data for ', num2str(N_points), '[V+H] observation points (Lebedev degree ', num2str(lebedevDegrees(degIdx)), ')']);
end

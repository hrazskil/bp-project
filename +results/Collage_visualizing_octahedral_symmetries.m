clc; clear; close all;

% Define different resolutions for Lebedev grid (example values)
N_values = [6, 14, 26];  % Low, Medium, High point counts
titles = {'Low Resolution', 'Medium Resolution', 'High Resolution'};
colors = {'r', 'g', 'b'};

% Octahedral symmetry transformations (24 elements)
octahedral_rotations = {
    eye(3), rotx(90), rotx(180), rotx(270), ...
    roty(90), roty(180), roty(270), ...
    rotz(90), rotz(180), rotz(270), ...
    rotx(90) * roty(90), rotx(90) * roty(-90), ...
    rotx(90) * rotz(90), rotx(-90) * rotz(90), ...
    rotx(90) * rotz(-90), rotx(-90) * rotz(-90), ...
    roty(45) * rotz(90), roty(-45) * rotz(90), ...
    roty(45) * rotz(-90), roty(-45) * rotz(-90), ...
    -eye(3) % Inversion symmetry
};
rotation_titles = {
    'Identity', 'Rotate 90° X', 'Rotate 180° X', 'Rotate 270° X', ...
    'Rotate 90° Y', 'Rotate 180° Y', 'Rotate 270° Y', ...
    'Rotate 90° Z', 'Rotate 180° Z', 'Rotate 270° Z', ...
    'Rotate 90° X and Y', 'Rotate 90° X and -Y', 'Rotate 90° X and Z', ...
    'Rotate -90° X and Z', 'Rotate 90° X and -Z', 'Rotate -90° X and -Z', ...
    'Rotate 45° Y and Z', 'Rotate -45° Y and Z', 'Rotate 45° Y and -Z', 'Rotate -45° Y and -Z', ...
    'Inversion' % Inversion symmetry
};
num_symmetry_steps = length(octahedral_rotations);

% Select a reference point index to track (e.g., first point in each grid)
ref_index = 1;

% Axis limits for better visibility (adjust as needed)
axis_limit = 1.5;

% Loop through each resolution (page)
for i = 1:3
    % Create a new figure for each page (resolution)
    figure;
    hold on;
    axis off;
    set(gcf, 'Position', [100, 100, 1200, 800]);
    
    % Loop through symmetry transformations
    for step = 1:num_symmetry_steps
        % Generate Lebedev points for the current resolution
        [xyz, ~] = utilities.getLebedevSphere(N_values(i)); % Replace with real Lebedev function
        
        % Get the symmetry transformation matrix
        symmetry_transform = octahedral_rotations{step};
        
        % Apply octahedral symmetry transformation
        xyz_rot = xyz * symmetry_transform';
        
        % Reference points
        ref_point_fixed = xyz(ref_index, :); % Fixed reference point (original position)
        ref_point_moving = xyz_rot(ref_index, :); % Moving reference point after symmetry transformation
        
        % Create a new subplot for this symmetry transformation
        subplot(6, 4, step); % Arrange as 6 rows and 4 columns
        hold on; grid on; axis equal;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('%s: %s', titles{i}, rotation_titles{step}));
        
        % Set consistent axis limits for visibility
        axis([-axis_limit, axis_limit, -axis_limit, axis_limit, -axis_limit, axis_limit]);
        
        % Plot rotated grid points
        scatter3(xyz_rot(:,1), xyz_rot(:,2), xyz_rot(:,3), 50, colors{i}, 'filled');
        
        % Plot reference points
        scatter3(ref_point_fixed(1), ref_point_fixed(2), ref_point_fixed(3), 100, 'w', 'filled'); % Fixed (White)
        scatter3(ref_point_moving(1), ref_point_moving(2), ref_point_moving(3), 100, 'y', 'filled'); % Moving (Yellow)
        
        % Display current grid size (number of points)
        text(0, 0, axis_limit * 0.9, ...
            sprintf('Grid Points: %d', size(xyz, 1)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'k');
        
        hold off;
    end
    
    % Adjust figure layout (optional)
    sgtitle(['Symmetry Transformations for ' titles{i}]); % Title for the entire page
end
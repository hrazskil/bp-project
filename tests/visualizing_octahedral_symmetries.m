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

% Animation loop: Apply each symmetry transformation
for i = 1:3
    % Create a new figure for each resolution
    figure;
    hold on; grid on; axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(titles{i});
    
    % Set consistent axis limits for visibility
    axis([-axis_limit, axis_limit, -axis_limit, axis_limit, -axis_limit, axis_limit]);
    
    % Set a consistent view angle (adjust if needed)
    view(3); % 3D view
    rotate3d on; % Enable rotation with mouse
    
    % Loop through symmetry transformations
    for step = 1:num_symmetry_steps
        symmetry_transform = octahedral_rotations{step}; % Get transformation matrix
        
        cla; % Clear frame for smooth animation
        
        % Generate Lebedev points for the current resolution
        [xyz, ~] = utilities.getLebedevSphere(N_values(i)); % Replace with real Lebedev function
        
        % Apply octahedral symmetry transformation
        xyz_rot = xyz * symmetry_transform';
        
        % Reference points
        ref_point_fixed = xyz(ref_index, :); % Fixed reference point (original position)
        ref_point_moving = xyz_rot(ref_index, :); % Moving reference point after symmetry transformation
        
        % Plot rotated grid points
        scatter3(xyz_rot(:,1), xyz_rot(:,2), xyz_rot(:,3), 50, colors{i}, 'filled');
        
        % Plot reference points
        scatter3(ref_point_fixed(1), ref_point_fixed(2), ref_point_fixed(3), 100, 'w', 'filled'); % Fixed (White)
        scatter3(ref_point_moving(1), ref_point_moving(2), ref_point_moving(3), 100, 'y', 'filled'); % Moving (Yellow)
        
        % Display current rotation and grid size
        text(0, 0, axis_limit * 0.9, ...
            sprintf('Rotation: %s\nGrid Points: %d', rotation_titles{step}, size(xyz, 1)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'k');
        
        % Pause to observe each symmetry transformation
        pause(1); 
    end
end

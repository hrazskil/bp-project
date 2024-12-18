function [r, theta, phi] = cart2sphCoor(x, y, z)
    % CART2SPHCOOR Converts Cartesian coordinates to spherical coordinates
    % Input:
    %  x     - x-coordinate in Cartesian system
    %  y     - y-coordinate in Cartesian system
    %  z     - z-coordinate in Cartesian system
    % Output:
    %  r     - radius (distance from origin)
    %  theta - polar angle (angle from the z-axis)
    %  phi   - azimuthal angle (angle in the x-y plane from the x-axis)

    % Calculate the radius
    r = sqrt(x.^2 + y.^2 + z.^2); % Radial distance from the origin

    % Calculate the polar angle theta
    theta = angle(z + 1i * sqrt(x.^2 + y.^2));
    theta(isnan(theta)) = 0; % Set NaN values to 0 for theta

    % Calculate the azimuthal angle phi
    phi = angle(x + 1i * y);
    phi(isnan(phi)) = 0; % Set NaN values to 0 for phi
end
function [x,y,z] = sph2cartCoor(r,theta,phi)

% SPH2CARTCOOR Converts spherical coordinates to cartesian coordinates
    % Input:
    %  r     - radius (distance from origin)
    %  theta - polar angle (angle from the z-axis)
    %  phi   - azimuthal angle (angle in the x-y plane from the x-axis)
    % Output:
    %  x     - x-coordinate in Cartesian system
    %  y     - y-coordinate in Cartesian system
    %  z     - z-coordinate in Cartesian system
    
    % Calculate the x-component of the coordinates
    x       = r.*cos(phi).*sin(theta);
    % Calculate the y-component of the coordinates
    y       = r.*sin(phi).*sin(theta);
    % Calculate the z-component of the coordinates
    z       = r.*cos(theta);
end


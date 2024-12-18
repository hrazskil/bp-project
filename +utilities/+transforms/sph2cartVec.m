function [F_x, F_y, F_z] = sph2cartVec(theta, phi, F_r, F_theta, F_phi)
    % SPH2CARTVEC Converts spherical forces to Cartesian coordinates.
    % Inputs:
    %  theta   - Azimuthal angle (angle from the z-axis)
    %  phi     - Polar angle (angle in the x-y plane from the x-axis)
    %  F_r     - Radial vector component
    %  F_theta - Azimuthal vector component
    %  F_phi   - Polar vector component
    % Outputs:
    %  F_x,F_y,F_z - Components of the vector in Cartesian coordinates

    % Repeated calculations of sine and cosine for angles
    sthe = sin(theta); % Sine of theta
    cthe = cos(theta); % Cosine of theta
    sphi = sin(phi);   % Sine of phi
    cphi = cos(phi);   % Cosine of phi

    % Calculate the x-component of the vector
    F_x = F_r .* sthe .* cphi + F_theta .* cthe .* cphi - F_phi .* sphi; 

    % Calculate the y-component of the vector
    F_y = F_r .* sthe .* sphi + F_theta .* cthe .* sphi + F_phi .* cphi;

    % Calculate the z-component of the vector
    F_z = F_r .* cthe - F_theta .* sthe; 
end




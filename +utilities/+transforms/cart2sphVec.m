function [F_r, F_theta, F_phi] = cart2sphVec(x, y, z, F_x, F_y, F_z)
    % CART2SPHVEC Converts Cartesian coordinates and vector components to spherical coordinates.
    % Inputs:
    %  x,y,z        - Cartesian coordinates
    %  F_x,F_y,F_z  - Components of the vector in Cartesian coordinates
    % Outputs:
    %  F_r      - Radial component of the vector
    %  F_theta  - Polar (theta) component of the vector
    %  F_phi    - Azimuthal (phi) component of the vector
  
    % Repeated calculations for radius and xy-plane distance
    r = sqrt(x.^2 + y.^2 + z.^2);   % Radial distance from the origin
    rho = sqrt(x.^2 + y.^2);         % Distance in the xy-plane
    phi = angle(x + 1i*y); % spherical angle phi

    % Calculate the radial component of the vector
    F_r = (x.*F_x + y.*F_y + z.*F_z) ./ r;
    
    % Calculate the Polar component of the vector
    % F_theta = (F_x.*z.*x + F_y.*z.*y - F_z.*rho.^2) ./ (r.*rho);
    F_theta = (F_x.*z.*cos(phi) + F_y.*z.*sin(phi) - F_z.*rho)./r;

    % Calculate the azimuthal component of the vector
    % F_phi = (-F_x.*y + F_y.*x) ./ rho;
    F_phi = (-F_x.*sin(phi) + F_y.*cos(phi)) ;

end
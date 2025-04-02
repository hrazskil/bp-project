function [powerDensity] = powerDensityFar(eF)
% POWERDENSITYFAR 
% Computes the far-field power density from the electric field vectors.
%
% Inputs:
%   eF        - Nx3 matrix of electric field vectors at observed points (V/m)
%
% Outputs:
%   powerDensity - Nx1 vector of power densities (W/m²) at the observed points

    % Retrieve constants
    const = utilities.constants.giveConstants; 
    Z0 = const.Z0;  % Free-space impedance (Ω)

    %% Step 1: Compute Power Density (Scalar)
    % Sum the squared magnitudes of the electric field components at each point
    % and scale by 1 / (2 * Z0)
    powerDensity = (1 / (2 * Z0)) * sum(abs(eF).^2, 2);  % Nx1 vector

end

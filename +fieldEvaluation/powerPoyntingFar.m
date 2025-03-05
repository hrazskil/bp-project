function [powerPoyntFar] = powerPoyntingFar(eF)
% POWERPOYNTINGFAR 
% Computes the far-field power density from the electromagnetic field.
%   Inputs:
%       eF - Nx3 matrix of electric field vectors at observed points.
%   Output:
%       pow - Power density at observed points (W/m^2).

    % Retrieve constants
    const = utilities.constants.giveConstants; 
    
    % Calculate the power density
    powerPoyntFar = (1 / (2 * const.Z0)) .* sum(abs(eF).^2, 2);
end

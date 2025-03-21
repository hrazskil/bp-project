function [integral] = powerQuadratureFar(Nleb, dip, f)
    % powerQuadratureFar computes the total radiated power using Lebedev quadrature
    % Inputs:
    %   Nleb - Number of Lebedev points used for quadrature
    %   dip  - Structure containing dipole parameters
    %   f    - Frequency of operation
    % Output:
    %   integral - The total radiated power computed using quadrature

    % Get Lebedev points and weights for numerical integration
    [points, weights, ~] = utilities.getLebedevSphere(Nleb);

    % Get physical constants
    construct = utilities.constants.giveConstants;
    
    % Compute the far-field electric field at Lebedev points
    fF = fieldEvaluation.farField(points, dip, f);

    % Compute the radiated power integral using Lebedev quadrature
    integral = sum(sum(fF .* conj(fF), 2) .* weights) / (2 * construct.Z0);
end


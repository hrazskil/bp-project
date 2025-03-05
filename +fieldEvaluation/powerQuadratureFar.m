function [integral] = powerQuadratureFar(Nleb,dip,f)
    % powerQuadrature computes the integral using Lebedev quadrature.
    % Inputs:
    %   Nleb - Number of Lebedev points
    %   eF   - Input vector for the electric field at observation points
    %   mf   - Input vector for the magnetic field at observation points
    % Output:
    %   integral - The integrated power radiated

    % Get Lebedev points and weights for the specified number of points
    [points, weights, ~] = utilities.getLebedevSphere(Nleb);

    % get constants
    construct   = utilities.constants.giveConstants;
    
    fF = fieldEvaluation.farField(points,dip,f);

    % Compute the integral
    integral = sum(sum(fF.*conj(fF),2).* weights)/(2*construct.Z0);
end


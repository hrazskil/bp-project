function [integral] = powerQuadrature(Nleb, dip, f)
    % powerQuadrature computes the integral using Lebedev quadrature.
    % Inputs:
    %   Nleb - Number of Lebedev points
    %   eF   - Input vector for the electric field at observation points
    %   mf   - Input vector for the magnetic field at observation points
    % Output:
    %   integral - The integrated power radiated
    %% quadrature
    
    % Get Lebedev points and weights for the specified number of points
    [rObserved, weights, ~] = utilities.getLebedevSphere(Nleb);

    [eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f);
    [mF] = fieldEvaluation.magFieldM2(rObserved,dip,f);

    % Calculate the temporary value
    tmp = real(sum((utilities.rowCross(eF, conj(mF)) .* rObserved), 2));

    % Compute the final integral by summing the product of tmp and weights
    integral = 0.5 * sum(tmp .* weights);
end


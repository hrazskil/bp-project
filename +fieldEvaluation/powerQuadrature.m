function [integral] = powerQuadrature(Nleb, eF, mF)
    % powerQuadrature computes the integral using Lebedev quadrature.
    % Inputs:
    %   Nleb - Number of Lebedev points
    %   eF   - Input vector for the electric field at observation points
    %   mf   - Input vector for the magnetic field at observation points
    % Output:
    %   integral - The integrated power radiated

    % Get Lebedev points and weights for the specified number of points
    [points, weights, ~] = utilities.getLebedevSphere(Nleb);

    % Calculate the temporary value
    tmp = 0.5 * real(sum((utilities.rowCross(eF, conj(mF)) .* points), 2));

    % Compute the final integral by summing the product of tmp and weights
    integral = sum(tmp .* weights);
end


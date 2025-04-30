function error = normObjectiveFunction(ampGuess, dip, f, points, weights, fF_ref, totalPower_ref)
    % normObjectiveFunction computes the normalized error between estimated and reference far-field patterns.
    % Inputs:
    %   ampGuess       - Current estimate of dipole amplitudes (complex vector)
    %   dip            - Dipole structure with positions and directions
    %   f              - Frequency of operation
    %   points         - Number of Lebedev quadrature points
    %   weights        - Lebedev quadrature weights for numerical integration
    %   fF_ref         - Reference far-field pattern
    %   totalPower_ref - Reference total radiated power
    % Output:
    %   error          - Normalized squared difference (percentage of total radiated power)
    
    % Physical constants
    construct = utilities.constants.giveConstants();

    % Update dipole amplitudes with current guess
    dipOpt = dip;
    dipOpt.complAmpl = ampGuess;  % Already in complex format

    % Compute far-field pattern for estimated dipoles
    fF = fieldEvaluation.farFieldM2(points, dipOpt, f);

    totalPower = sum(weights .* sum(abs(fF).^2, 2)) / (2 * construct.Z0);

    % Compute squared difference and integrate over sphere
    error = sum(weights .* sum(abs(abs(fF)/sqrt(totalPower) - abs(fF_ref)/sqrt(totalPower_ref)).^2, 2)) / (2 * construct.Z0);

    % error = sum(weights .* sum(abs(fF/sqrt(totalPower) - fF_ref/sqrt(totalPower_ref)).^2, 2)) / (2 * construct.Z0);
end

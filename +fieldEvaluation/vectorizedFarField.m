function [fF] = vectorizedFarField(rObserved, dip, f)
    % Vectorized function to compute the far-field radiation pattern from dipoles
    nDip = size(dip.pos, 1); % Number of dipoles
    nObs = size(rObserved, 1); % Number of observed points

    % Constants
    construct = utilities.constants.giveConstants; % Retrieve physical constants
    omega = 2 * pi * f; % Angular frequency
    k = omega / construct.c0; % Wave number (k = omega/c)

    % Normalize observed positions
    r0 = rObserved ./ utilities.rowNorm(rObserved); % Normalize each observed position
    r0 = reshape(r0', 1, 3, nObs); % Reshape for broadcasting

    % Compute the multipole moment
    Mp = reshape(dip.complAmpl, 1, nDip) .* reshape(dip.dir', 3, 1, nDip); % Combine amplitude and direction

    % Compute the far-field pattern
    expTerm = exp(1i * k * sum(r0 .* reshape(dip.pos', 3, 1, nDip), 2)); % Exponential term
    crossTerm = cross(-r0, cross(r0, Mp, 2), 2); % Cross product term

    fF = construct.Z0 * construct.c0 * k^2 * expTerm ./ (4 * pi) .* crossTerm; % Final computation

    % Sum over dipoles and reshape
    fF = reshape(sum(fF, 3), 3, []).'; % Reshape to Mx3 matrix
end

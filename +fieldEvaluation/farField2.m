function [fF] = farField2(rObserved, dip, f)
    % FARFIELD Computes the far-field radiation pattern from dipoles.
    % Inputs:
    %   rObserved - Nx3 matrix of observed positions (x, y, z)
    %   dip - structure containing dipole properties (positions, amplitudes, directions)
    %   f - frequency of the source
    % Output:
    %   fF - far-field radiation pattern (Nx3 matrix)

    % Number of dipoles and observation points
    nDip = size(dip.pos, 1);
    nObs = size(rObserved, 1);

    % Retrieve physical constants
    constants = utilities.constants.giveConstants();
    omega = 2 * pi * f;                 % Angular frequency
    k = omega / constants.c0;           % Wave number (k = omega / c)
    Z0 = constants.Z0;                  % Impedance of free space

    % Normalize observed positions
    rNorm = utilities.rowNorm(rObserved); % Compute norms of each observed position
    r0 = rObserved ./ rNorm;           % Normalize each observed position
    r0 = repelem(r0, nDip, 1);         % Repeat for each dipole (Nx3 repeated nDip times)

    % Precompute repeated terms
    dipPosRep = repelem(dip.pos, nObs, 1);          % Repeat dipole positions for all observation points
    expTerm = exp(1i * k * sum(r0 .* dipPosRep, 2)); % Exponential term: exp(ik * (r0 · dip.pos))

    % Compute the multipole moment for each dipole
    Mp = repelem(dip.complAmpl, nObs, 1) .* repelem(dip.dir, nObs, 1);

    % Calculate the far-field components
    crossMpR0 = cross(-r0, cross(r0, Mp, 2), 2);    % Double cross-product
    fF = (Z0 * constants.c0 * k^2 / (4 * pi)) .* expTerm .* crossMpR0;

    % Sum over all dipoles
    fF = sum(reshape(fF, nObs, nDip, 3), 2);       % Reshape and sum over dipoles
    fF = squeeze(fF);                              % Remove singleton dimension (Mx3 matrix)
end

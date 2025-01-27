function [powerRad] = powerRadiated(f, dip)
% powerRadiated Calculates the total power radiated by the dipoles.
% Inputs:
%   f - Frequency of the source.
%   dip - Structure containing the dipole properties:
%       dip.complAmpl - Vector of complex amplitudes (dipole moments).
% Outputs:
%   powerRad - Total power radiated by all dipoles.

const = utilities.constants.giveConstants;
omega = 2*pi*f;               % Angular frequency
k = omega / const.c0;     % Wave number

% Total radiated power for all dipoles
powerRad = (const.c0^2 * const.Z0 * k^4 / (12 * pi)) * ...
            sum(abs(dip.complAmpl).^2, 1); 
% Sum of squared magnitudes of dipole moments
end

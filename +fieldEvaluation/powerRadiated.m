function [powerRad] = powerRadiated(f,dip)

% powerRad = calculated power radiated by the dipoles whose dipole
% moment magnitudes are submited via the ComplAmpl vector
% f = frequency
% ComplAmpl = vector containing the magnitudes of the dipole moments of individual dipoles

construct   = utilities.constants.giveConstants;
omega       = 2*pi*f;
k           = omega/construct.c0;

powerRad    = (construct.c0^2*construct.Z0*k^4/(12*pi))* ...
            sum(abs(dip.complAmpl).^2,1);
end


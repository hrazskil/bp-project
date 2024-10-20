function [powerRad] = powerRadiated(f,complAmpl)

%POWERRADIATED calculation of the full power radiated by the dipoles
% f = frequency

construct   = utilities.constants.giveConstants;
omega       = 2*pi*f;
k           = omega/construct.c0;

powerRad    = (construct.c0^2*construct.Z0*k^4/(12*pi))* ...
            sum(abs(complAmpl).^2,1);
end


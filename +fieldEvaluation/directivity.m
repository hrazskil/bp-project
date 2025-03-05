function [directivity] = directivity(Imax,totalPower)
% Directivity calculation
directivity = (4 * pi * Imax) / totalPower; % Directivity formula
end


function [powerPoynt] = powerPoynting(eF,mF)
% POWERPOYNTING
%   Computes the power density of the electromagnetic field at observation 
% points.
%   Inputs:
%       eF - Nx3 matrix of electric field vectors at observed points.
%       mF - Nx3 matrix of magnetic field vectors at observed points.
%   Output:
%       powerPoynt - Nx3 matrix of the real part of the Poynting vector.

    % Compute the real part of the Poynting vector: S = 0.5 * Re(E × H*)
    powerPoynt = real(utilities.rowCross(eF, conj(mF))) / 2;
end


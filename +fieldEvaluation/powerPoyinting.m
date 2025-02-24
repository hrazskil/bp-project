function [powerPoyint] = powerPoyinting(eF,mF)
% POWERPOYNTING 
% Computes the power density of the electromagnetic field at observation 
% points.
%   Inputs:
%       eF - Nx3 matrix of electric field vectors at observed points.
%       mF - Nx3 matrix of magnetic field vectors at observed points.
%   Output:
%       pow - Power density at observed points (W/m^2).

    % Calculate the power density
    powerPoyint = real(utilities.rowCross(eF, conj(mF)))/2;
end


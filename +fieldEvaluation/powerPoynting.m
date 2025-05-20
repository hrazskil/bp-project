function [powerPoynt] = powerPoynting(eF, mF)
% powerPoynting Computes time-averaged Poynting vector from E and H fields.
%
% INPUTS:
%   eF (nObs × 3) - Electric field vectors (complex).
%   mF (nObs × 3) - Magnetic field vectors (complex).
%
% OUTPUT:
%   powerPoynt (nObs × 3) - Real part of the time-averaged Poynting vector.
%
% FORMULATION:
%   S = 0.5 * Re(E × H*)
%
% ------------------------------------------------------------------------

% Step 1: Compute the real part of E × H*
powerPoynt = (1/2)*real(utilities.rowCross(eF, conj(mF))) / 2;

end

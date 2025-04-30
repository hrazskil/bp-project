function [integral] = powerQuadrature(Nleb, dip, f)
% powerQuadrature Computes total radiated power using Poynting vector integration.
%
% INPUTS:
%   Nleb (scalar) - Number of Lebedev quadrature points.
%   dip  (struct) - Dipole configuration (positions, amplitudes, directions).
%   f    (scalar) - Frequency in Hz.
%
% OUTPUT:
%   integral (scalar) - Estimated total radiated power.
%
% FORMULATION:
%   Integral over a unit sphere of Re(E × H*) • r_hat, weighted by quadrature.
%
% ------------------------------------------------------------------------

% Step 1: Get Lebedev points and weights
[rObserved, weights, ~] = utilities.getLebedevSphere(Nleb);

% Step 2: Evaluate electric and magnetic fields at these points
eF = fieldEvaluation.eleFieldM2(rObserved, dip, f);
mF = fieldEvaluation.magFieldM2(rObserved, dip, f);

% Step 3: Compute Poynting vector component in radial direction
S_radial = real(sum(utilities.rowCross(eF, conj(mF)) .* rObserved, 2));

% Step 4: Integrate using weighted sum
integral = 0.5 * sum(S_radial .* weights);

end

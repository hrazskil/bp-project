function [propagationFactor] = computePropagationFactor(rObserved, k)
% COMPUTEELECTRICFIELD Converts far-field radiation pattern F to electric field E.
% Inputs:
%   rObserved - Nx3 matrix of observed positions (x, y, z).
%   k - Wave number (2*pi*f / c).
% Output:
%   eF - Electric field at observed points (Nx3 matrix).

% Compute distances (r) 
r = utilities.rowNorm(rObserved); %
% Apply the propagation factor to far-field F
propagationFactor = exp(-1i * k * r) ./ r; % Nx1 vector
end
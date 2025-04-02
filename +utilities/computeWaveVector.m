function kVec = computeWaveVector(rObserved, dipPos, freq)
% COMPUTEWAVEVECTOR 
% Computes the wave vector for far-field calculations.
%
% Inputs:
%   rObserved - Nx3 matrix of observation points (m)
%   dipPos    - Mx3 matrix of dipole positions (m)
%   freq      - Scalar frequency (Hz)
%
% Outputs:
%   kVec - Nx3 matrix of wave vectors (1/m)

    %% Step 1: Define Constants
    construct = utilities.constants.giveConstants; % Retrieve fundamental constants
    c0 = construct.c0; % Speed of light in vacuum (m/s)
    
    omega = 2 * pi * freq; % Angular frequency (rad/s)
    k_mag = omega / c0;    % Wavenumber k (1/m)

    %% Step 2: Compute Relative Position Vectors
    rVec = rObserved - dipPos;   % Nx3 matrix

    %% Step 3: Normalize to Get Unit Direction Vectors
    k_hat = rVec ./ utilities.rowNorm(rVec); % Normalize each row

    %% Step 4: Compute Wave Vector
    kVec = k_mag * k_hat;  % Multiply by wavenumber
    
end
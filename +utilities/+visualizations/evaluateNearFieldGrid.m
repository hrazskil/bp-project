function [E, H, S, grid] = evaluateNearFieldGrid(dip, freq, gridSpec)
% evaluateNearFieldGrid computes the near-field electric and magnetic fields,
% as well as the time-averaged Poynting vector, on a user-defined Cartesian grid.
%
% INPUTS:
%   dip      (struct) - Structure containing dipole model:
%       dip.pos       (nDip × 3)  - Cartesian positions of the dipoles.
%       dip.dir       (nDip × 3)  - Unit vectors representing dipole orientations.
%       dip.complAmpl (nDip × 1)  - Complex dipole moments (amplitudes).
%
%   freq     (scalar) - Operating frequency in Hz.
%
%   gridSpec (struct) - Structure defining the evaluation region:
%       gridSpec.x (1 × Nx) - Vector of x-coordinates.
%       gridSpec.y (1 × Ny) - Vector of y-coordinates.
%       gridSpec.z (1 × Nz) - Vector of z-coordinates.
%
% OUTPUTS:
%   E    (Ny × Nx × Nz × 3) - Complex electric field vectors at each grid point.
%   H    (Ny × Nx × Nz × 3) - Complex magnetic field vectors at each grid point.
%   S    (Ny × Nx × Nz)     - Magnitude of time-averaged Poynting vector.
%   grid (struct)           - Struct containing meshgrid coordinate arrays:
%       grid.X, grid.Y, grid.Z (Ny × Nx × Nz)
%
% FORMULATION:
% The function evaluates near-field E and H using the full analytical dipole expressions,
% then computes the time-averaged Poynting vector via:
%
%   S(r) = (1/2) * Re{ E(r) × conj(H(r)) }
%
% ------------------------------------------------------------------------

%% Step 2: Generate Cartesian Grid
[X, Y, Z] = meshgrid(gridSpec.x, gridSpec.y, gridSpec.z);
gridPts = [X(:), Y(:), Z(:)]; % Flattened list of evaluation points

%% Step 3: Evaluate Near-Field Vectors
Eflat = fieldEvaluation.eleFieldM2(gridPts, dip, freq); % (Npts × 3)
Hflat = fieldEvaluation.magFieldM2(gridPts, dip, freq); % (Npts × 3)

%% Step 4: Reshape Field Data into 4D Arrays
Ny = length(gridSpec.y);
Nx = length(gridSpec.x);
Nz = length(gridSpec.z);
E = reshape(Eflat, [Ny, Nx, Nz, 3]);
H = reshape(Hflat, [Ny, Nx, Nz, 3]);

%% Step 5: Compute Poynting Vector Magnitude
%   S = (1/2) * Re{E × conj(H)}
Sx = 0.5 * real(E(:,:,:,2).*conj(H(:,:,:,3)) - E(:,:,:,3).*conj(H(:,:,:,2)));
Sy = 0.5 * real(E(:,:,:,3).*conj(H(:,:,:,1)) - E(:,:,:,1).*conj(H(:,:,:,3)));
Sz = 0.5 * real(E(:,:,:,1).*conj(H(:,:,:,2)) - E(:,:,:,2).*conj(H(:,:,:,1)));
S = sqrt(Sx.^2 + Sy.^2 + Sz.^2); % Magnitude of total Poynting vector

%% Step 6: Return Grid Struct
grid.X = X;
grid.Y = Y;
grid.Z = Z;

end
function [eM] = magFieldM2(rObserved,dip,f)
% MAGFIELDM2 Computes the magnetic field at observation points due to dipoles.
% 
% Inputs:
%   rObserved - Nx3 matrix of observation points in space.
%   dip       - Struct with dipole properties:
%               dip.pos (Mx3)      - Positions of dipoles.
%               dip.dir (Mx3)      - Unit direction vectors of dipoles.
%               dip.complAmpl (Mx1) - Complex amplitudes of dipoles.
%   f         - Frequency of operation (Hz).
% 
% Output:
%   eM - Nx3 matrix representing the magnetic field at observation points.
%
% Notes:
%   - Uses vectorized computations for efficiency.


    % Number of dipoles and observation points
    nDip    = size(dip.pos,1);
    nObs    = size(rObserved,1);
    
     % Physical constants
    construct = utilities.constants.giveConstants;  % Load constants
    omega   = 2*pi*f;                               % Angular frequency
    k       = omega/construct.c0;                   % Wavenumber (k = ω/c)
    
    % Compute displacement vectors from dipoles to observation points
    MRVec   = repmat(rObserved,nDip,1) - repelem(dip.pos,nObs,1); % (MxN)x3
    MR      = utilities.rowNorm(MRVec); % Compute the Euclidean distance
    MdirR   = MRVec./MR;                % Normalize displacement vectors
    
    % Compute dipole moments at observation points
    Mp      = repelem(dip.complAmpl,nObs,1) .* (repelem(dip.dir,nObs,1));
    
    % Compute phase term (1i * k * R)
    IMR     = 1i*k*MR;
    
    % Compute the magnetic field using the dipole radiation formula
    eM     = construct.c0*k^2*exp(-IMR) ./ (4*pi*MR) .* (...
        utilities.rowCross(MdirR,Mp).*(1 ./ (IMR) + 1) );
    
    % Sum contributions from all dipoles and reshape output
    eM     = reshape(sum(reshape(eM', [], nDip), 2),3,[])';
end
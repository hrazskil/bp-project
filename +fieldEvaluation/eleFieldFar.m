function [eF] = eleFieldFar(rObserved, dip, f)
% ELEFIELDFAR aproximates the electric field in the far field.

% Initialization
nDip = size(dip.pos, 1);
nObs = size(rObserved, 1);

% Constants
construct = utilities.constants.giveConstants; 
omega = 2 * pi * f; 
k = omega / construct.c0; 

% Preallocate and calculate differences (it helps)
dipPosRep = repelem(dip.pos, nObs, 1); 
obsPosRep = repmat(rObserved, nDip, 1);
MRVec = obsPosRep - dipPosRep; 
MR = utilities.rowNorm(MRVec);

% Normalize directional vectors
MdirR = MRVec ./ MR;

% Compute dipole moment
Mp = repelem(dip.complAmpl, nObs, 1) .* repelem(dip.dir, nObs, 1);

% Compute the electric field
eF = construct.Z0 * construct.c0 * k^2 .* exp(-1i * k * MR) ./ (4 * pi * MR) ...
    .* (utilities.rowCross(-MdirR, utilities.rowCross(MdirR, Mp)));

% Sum contributions and reshape (more clear)
eF = sum(reshape(eF, nObs, nDip, 3), 2);
eF = reshape(eF, nObs, 3);

end
function [eF] = eleFieldM2(rObserved,dip,f)
%   eleField Summary of this function goes here
%   Detailed explanation goes here
% Initialization
nDip    = size(dip.pos,1);
nObs    = size(rObserved,1);

% Constants
construct = utilities.constants.giveConstants;
omega   = 2*pi*f;
k       = omega/construct.c0;


MRVec   = repmat(rObserved,nDip,1) - repelem(dip.pos,nObs,1);
MR      = utilities.rowNorm(MRVec);
MdirR   = MRVec./MR;

Mp      = repelem(dip.complAmpl,nObs,1) .* (repelem(dip.dir,nObs,1));


eF     = construct.Z0*construct.c0*k^2.*exp(-1i*k*MR)./(4*pi*MR).*...
    (...
    utilities.rowCross(-(MdirR),utilities.rowCross((MdirR),Mp)) +...
        (...
        (3*(MdirR).*(utilities.rowDot(MdirR,Mp))-Mp) .* ...
        (1./(k^2*MR.^2)+(1i./(k*MR)))...
        )...
    );

% Sum contributions from all dipoles
eF     = reshape(sum(reshape(eF', [], nDip), 2),3,[])';
end
    

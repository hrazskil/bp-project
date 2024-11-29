function [eF] = eleFieldM2(rObserved,dip,f,complAmpl)
%   eleField Summary of this function goes here
%   Detailed explanation goes here
%inicialization
nDip    = size(dip.pos,1);
nObs    = size(rObserved,1);
eM2     = nan(nObs,3);

%constants
construct = utilities.constants.giveConstants;
omega   = 2*pi*f;
k       = omega/construct.c0;


MRVec   = repmat(rObserved,nDip,1) - repelem(dip.pos,nObs,1);
MR      = utilities.rowNorm(MRVec);
MdirR   = MRVec./MR;

Mp      = repelem(complAmpl,nObs,1) .* (repelem(dip.dir,nObs,1));


eF     = construct.Z0*construct.c0*k^2.*exp(-1i*k*MR)./(4*pi*MR).*(...
    cross(-(MdirR),cross((MdirR),Mp,2),2) +...
    ( (3*(MdirR).*(dot(MdirR,Mp,2))-Mp).* (1./(k*MR.^2)+(1i./(MR))) ) );

% Hell of a sum
eF     = reshape(sum(reshape(eM2', [], nDip), 2),3,[])';
end
    

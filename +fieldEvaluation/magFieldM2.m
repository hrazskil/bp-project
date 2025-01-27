function [eM] = magFieldM2(rObserved,dip,f)
%MAGFIELDM Summary of this function goes here
%inicialization
nDip    = size(dip.pos,1);
nObs    = size(rObserved,1);

%constants
construct = utilities.constants.giveConstants;
omega   = 2*pi*f;
k       = omega/construct.c0;


MRVec   = repmat(rObserved,nDip,1) - repelem(dip.pos,nObs,1);
MR      = utilities.rowNorm(MRVec);
MdirR   = MRVec./MR;

Mp      = repelem(dip.complAmpl,nObs,1) .* (repelem(dip.dir,nObs,1));

IMR     = 1i*k*MR;

eM     = construct.c0*k^2*exp(-IMR) ./ (4*pi*MR) .* (...
    utilities.rowCross(MdirR,Mp).*(1 ./ (IMR) + 1) );

% Hell of a sum
eM     = reshape(sum(reshape(eM', [], nDip), 2),3,[])';

end


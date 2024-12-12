function [fF] = farField(rObserved,dip,f)
%ELEFIELDFAR Summary of this function goes here
%   Detailed explanation goes here
%inicialization
nDip    = size(dip.pos,1);
nObs    = size(rObserved,1);

%constants
construct = utilities.constants.giveConstants;
omega   = 2*pi*f;
k       = omega/construct.c0;


r0 = rObserved./repmat(utilities.rowNorm(rObserved),[1,3]);
r0 = repmat(r0,[nDip,1]);

Mp      = repelem(dip.complAmpl,nObs,1) .* (repelem(dip.dir,nObs,1));

fF     = construct.Z0*construct.c0*k^2.*exp(1i*k*sum(r0 .* repelem(dip.pos,nObs,1),2))./(4*pi).*(...
    cross(-r0,cross(r0,Mp,2),2) );

% Hell of a sum
fF     = reshape(sum(reshape(fF', [], nDip), 2),3,[])';

end


function [eF] = eleField(rObserved,dip,f,complAmpl)
%   eleField Summary of this function goes here
%   Detailed explanation goes here

%inicialization
construct   = utilities.constants.giveConstants;
nObs        = size(rObserved,1);
eF          = nan(nObs,3);

%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;

    for iObs = 1:nObs

    Rvec    = repmat(rObserved(iObs,:),[size(dip.pos,1),1])-dip.pos;
    R       = repmat(utilities.rowNorm(Rvec),[1,3]);
    dirR    = Rvec./R;
    
    p       = repmat(complAmpl,[1,3]).*dip.dir;

    eF(iObs,:) = sum(construct.Z0*construct.c0*k^2* ...
    exp(-1i*k*R)./(4*pi*R) .*( cross(-(dirR),cross((dirR),p,2),2) + ...
    ( (3*(dirR).*(dot(dirR,p,2))-p) .* (1./(k^2*R.^2)+(1i./(k*R))) ) ),1);
    end
end
    

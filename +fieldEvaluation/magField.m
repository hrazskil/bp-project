function [mF] = magField(rObserved,dip,f,complAmpl)
%MAGFIELD Summary goes here
%inicialization
construct   = utilities.constants.giveConstants;
nObs        = size(rObserved,1);
mF          = nan(nObs,3);

%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;

    for iObs = 1:nObs
    Rvec    = repmat(rObserved(iObs,:),[size(dip.pos,1),1])-dip.pos;
    R       = repmat(utilities.rowNorm(Rvec),[1,3]);
    dirR    = Rvec./R;
    p       = repmat(complAmpl,[1,3]).*dip.dir;

    mF(iObs,:) = sum(exp(-1i*k*R)./(4*pi*R) .*(...
        cross(dirR,p,2).*(1-1*i./(k*R))),1);
    end
 mF = construct.c0*k^2*mF;
end


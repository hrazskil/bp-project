function [outputArg1] = eleFieldFar(rObserved,rSource,phase,dipoleMomentSize)
%ELEFIELDFAR Summary of this function goes here
%   Detailed explanation goes here
c0=299792458;
k=phase/c0;
R=rObserved-rSource;
nR=norm(R);
dirR=R/nR;
eleFieldFar=-dipoleMomentSize*c0*k^3*
outputArg2 = inputArg2;
end


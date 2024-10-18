function [eZ] = eleField(rObserved,rSource,dipoleMoment,phase)
%*eleField Summary of this function goes here
%   Detailed explanation goes here
c0=299792458;
tic
c0= 3.0e+08;
k=phase/c0;
R=rObserved-rSource;
nR=norm(R);
dirR=R./nR;
eZ=[0;0;1].*( ...
    c0*k^3*exp(-1i*k*nR)/(4*pi*k*nR)*( ...
                                        cross(-(dirR),cross((dirR),dipoleMoment)) + ...
                                        (3*(dirR)*(dot(dipoleMoment,dirR)) -dipoleMoment) * (1/(k^2*nR^2) + 1i*(1/(k*nR))) ...
                                     )  ...
            )
toc
end


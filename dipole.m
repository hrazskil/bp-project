clc
clear
dip.pos=[0;0;0];
dip.phase = pi/4;
dip.mom=[0;0;1];




rObserved=[0;0;1];

phase = dip.phase;
rSource = dip.pos;
dipoleMoment= dip.mom;
tic

%*eleFieldSpherical
% c0=3.0e+08;
% k=phase/c0;
% R=(rObserved-rSource);
% nR=norm(R);
% dirR=R./nR;
% [~,theta,r] = cart2sph(R(1,:),R(2,:),R(3,:));
% eZ= dipoleMoment.*(c0*k^3*( dirR.*(cos(theta)*2*( 1/(k^2*r^2) + 1i/(k*r)) ) + [0; 1; 0]*sin(theta)*( -1 + 1/(k^2*r^2) + 1i/(k*r) ) )*exp(-1i*k*r)/(4*pi*k*r))
% toc
%ELEFIELD Summary of this function goes here
%   Detailed explanation goes here

tic
c0= 3.0e+08;
k=phase/c0;
R=rObserved-rSource;
nR=norm(R);
dirR=R./nR;
e=[0;0;1].*( ...
    c0*k^3*exp(-1i*k*nR)/(4*pi*k*nR)*( ...
                                        cross(-(dirR),cross((dirR),dipoleMoment)) + ...
                                        (3*(dirR)*(dot(dipoleMoment,dirR)) -dipoleMoment) * (1/(k^2*nR^2) + 1i*(1/(k*nR))) ...
                                     )  ...
            )
toc
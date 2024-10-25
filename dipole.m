clc
clear
nDip=20;
dip.pos = rand(nDip,3);
dip.dir = rand(nDip,3);
complAmpl = rand(nDip,1)+  1i*rand(nDip,1);


f = 10000;

nObs=20000;
rObserved=rand(nObs,3);
%%
% tic
% [eZ] = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
% toc
%
 % better for ndip=nobs or ndip>nobs
tic
[eM] = fieldEvaluation.magField(rObserved,dip,f,complAmpl);
toc

%better for nobs>>ndip
%needs to be optimized
tic
[eM2] = fieldEvaluation.magFieldM(rObserved,dip,f,complAmpl); 
toc
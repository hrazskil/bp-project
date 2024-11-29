%% cleaning
clc
clear
%% inicialization for testing 

nDip=2;
dip.pos     = rand(nDip,3);
dip.dir     = rand(nDip,3);
complAmpl   = rand(nDip,1)+  1i*rand(nDip,1);


f           = 10000;

nObs        = 10;
rObserved   = rand(nObs,3);

    %% tests of functions
tic
[eF] = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
toc

eF.*conj(eF)

tic
const = utilities.constants.giveConstants;
pow = (1/(2*const.Z0))*(eF.*conj(eF))
toc
%  % better for ndip=nobs or ndip>nobs
% tic
% [mF] = fieldEvaluation.magField(rObserved,dip,f,complAmpl);
% toc
% 
% %better for nobs>>ndip
% %needs to be optimized
% tic
% [mFM] = fieldEvaluation.magFieldM(rObserved,dip,f,complAmpl); 
% toc
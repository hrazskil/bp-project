%% cleaning
clc
clear

%% inicialization for testing functions
nDip=1;
dip.pos     = [0.4,0.3,0.7];
%%
dip.dir     = [1,2,3];
%% norming of dip.dir
dip.dir = dip.dir./repmat(utilities.rowNorm(dip.dir),[1,3]);
%%
complAmpl   = 0.4+1i*0.6;
%%
f           = 1e4;

%% control function

%inicialization
construct   = utilities.constants.giveConstants;

%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;

%% tests of functions




% tic
% const   = utilities.constants.giveConstants;
% pow     = (1/(2*const.Z0))*(eF.*conj(eF))
% toc    
%% functioning functions
% 
rObserved   = [1,2,3];
rObserved = (1000/k)*rObserved
%%
tic
[eF]    = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
toc

%  % better for ndip=nobs or ndip>nobs
tic
[mF] = fieldEvaluation.magFieldM2(rObserved,dip,f,complAmpl);
toc

%
% %better for nobs>>ndip
% %needs to be optimized
% tic
% [mFM] = fieldEvaluation.magFieldM(rObserved,dip,f,complAmpl);
% toc
% %better for nobs>>ndip
% %needs to be optimized
% tic
% [mFM2] = fieldEvaluation.magFieldM2(rObserved,dip,f,complAmpl);
% toc

%% quadrature
Nleb = 302;

% quadrature points and weights
% To integrate a function F(x,y,z) over a unit sphere, call
% sum(F(points).*weigths)
[points, weigths, ~] = utilities.getLebedevSphere(Nleb);

rObserved = points*1e3/k;

tic
[eF] = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
toc

tic
[eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f,complAmpl);
toc

tic
[mF] = fieldEvaluation.magFieldM2(rObserved,dip,f,complAmpl);
toc

tmp = 0.5*real(sum((cross( eF, conj(mF),2).*points ),2));


% numerical
integral    = sum( tmp .* weigths);
integral    = integral*(1e3/k)^2

result      = fieldEvaluation.powerRadiated(f,complAmpl)
error       = result-integral

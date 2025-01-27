% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%% inicialization for testing functions
nDip=1;
dip.pos     = [0.4,0.3,0.7];
%%
dip.dir     = [1,2,3];
%% norming of dip.dir
dip.dir = dip.dir./repmat(utilities.rowNorm(dip.dir),[1,3]);
%%
dip.complAmpl   = 0.4+1i*0.6;

%% control function

%inicialization
construct   = utilities.constants.giveConstants;

%out for-loop variables
omega       = 2*pi*f0List;
k           = omega/construct.c0;

%% tests of functions




% tic
% const   = utilities.constants.giveConstants;
% pow     = (1/(2*const.Z0))*(eF.*conj(eF))
% toc    
%% functioning functions
% 
rObserved   = [1,2,3];
%%
% tic
% [eF] = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
% toc

tic
[eF] = fieldEvaluation.eleFieldM2(rObservedFar,dip,f0List,complAmpl);
toc

tic
[eF] = fieldEvaluation.eleFieldFar(rObservedFar,dip,f0List,complAmpl);
toc

%  % better for ndip=nobs or ndip>nobs
% tic
% [mF] = fieldEvaluation.magFieldM2(rObserved,dip,f,complAmpl);
% toc

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

rFar = 1e6/k;
rObserved = points*rFar;

tic
[eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f0List);
toc

tic
[mF] = fieldEvaluation.magFieldM2(rObserved,dip,f0List);
toc

tmp = 0.5*real(sum((cross( eF, conj(mF),2).*points ),2));


[fF] = fieldEvaluation.farField(rObserved,dip,f0List);



% numerical
integral    = sum( tmp .* weigths);
integral    = integral*rFar^2

result      = fieldEvaluation.powerRadiated(f0List,dip)
error       = result-integral

sum(sum(fF.*conj(fF),2).* weigths)/(2*construct.Z0)

%% Epoxid
% tmp = 0.5*real(sum((cross( eF, conj(mF),2).*points ),2));
% 
% 
% % numerical
% integral    = sum( tmp .* weigths);
% integral    = integral*(1e3/k)^2
% 
% result      = fieldEvaluation.powerRadiated(f,complAmpl)
% error       = result-integral

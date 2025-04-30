%% --- 1. Initialization ---
clc; clear; close all;

% Load dipole variables from file 
load('C:\Users\kilia\Plocha\gitHub\bp-project\tests\test_structure_2\DipoleArray.mat');
f = f0List; % Frequency of 1 GHz
% dip.dir = dip.dir'
% dip.pos = dip.pos'
% dip.complAmpl = dip.complAmpl'
%% inicialization for testing functions
% nDip=1;
% dip.pos     = [0.4,0.3,0.7];
% dip.dir     = [1,2,3];
% dip.dir = dip.dir./repmat(utilities.rowNorm(dip.dir),[1,3]);
% dip.complAmpl   = 0.4+1i*0.6;

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
[rObserved] = utilities.getLebedevSphere(302);



% tic
% [eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f0List);
% toc
% 
% tic
% [eF2] = fieldEvaluation.eleFieldM3(rObserved,dip,f0List);
% toc
% 
% sum(abs(eF2)-abs(eF),"all")

tic
[mF] = fieldEvaluation.farFieldM2(rObserved,dip,f);
toc

tic
[mF2] = fieldEvaluation.farFieldM2(rObserved,dip,f);
toc

sum(abs(mF2)-abs(mF),"all")

% tic
% [mF] = fieldEvaluation.magFieldM2(rObserved,dip,f);
% toc
% 
% tic
% [mF2] = fieldEvaluation.magFieldM3(rObserved,dip,f);
% toc
% 
% sum(abs(mF2)-abs(mF),"all")
%% quadrature
% 
% 
% 
% 
% Nleb = 302;
% 
% % quadrature points and weights
% % To integrate a function F(x,y,z) over a unit sphere, call
% % sum(F(points).*weigths)
% [points, weigths, ~] = utilities.getLebedevSphere(Nleb);
% 
% rFar = 1e6/k;
% rObserved = points*rFar;
% 
% tic
% [eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f0List);
% toc
% 
% tic
% [mF] = fieldEvaluation.magFieldM2(rObserved,dip,f0List);
% toc
% 
% tmp = 0.5*real(sum((cross( eF, conj(mF),2).*points ),2));
% 
% 
% [fF] = fieldEvaluation.farField(rObserved,dip,f0List);
% 
% 
% 
% % numerical
% integral    = sum( tmp .* weigths);
% integral    = integral*rFar^2
% 
% result      = fieldEvaluation.powerDensityFar()
% error       = result-integral
% 
% sum(sum(fF.*conj(fF),2).* weigths)/(2*construct.Z0)

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

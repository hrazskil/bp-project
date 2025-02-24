% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%%
% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\halfwaveDipole.mat');
f=f0List;

% dip.pos = [1, 1, 0.1;
%            0.5, 0.3, 0.1;
%            0.1, 0.2, 0.3];
% dip.dir = [-1, 0, 2;
%            0, 3, 2;
%            4, 2, 0];
% dip.dir = dip.dir./utilities.rowNorm(dip.dir);
% dip.complAmpl = [1;
%                  2;
%                  3];
% f=1e6;
%%


%% quadrature

Nleb = 5810;

[points, weigths, ~] = utilities.getLebedevSphere(Nleb);

%inicialization
construct   = utilities.constants.giveConstants;

%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;
rFar = 1e8/k;
rObserved = points*rFar;
%rObserved = [1, 2, 3]*rFar;
tic
[eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f);


[mF] = fieldEvaluation.magFieldM2(rObserved,dip,f);
toc

tmp = 0.5*real(sum((utilities.rowCross( eF, conj(mF)).*points ),2));

tic
[fF] = fieldEvaluation.farField(rObserved,dip,f);
propFactor = fieldEvaluation.computePropagationFactor(fF, rObserved, k);
eF2 = fF .* propFactor;
toc

%numerical
integral    = sum( tmp .* weigths);
integral    = integral*rFar^2
integralTest = fieldEvaluation.powerQuadrature(Nleb, eF, mF)*rFar^2
integral2   = sum(sum(fF.*conj(fF),2).* weigths)/(2*construct.Z0)

%sum of all errors
allErrorSumEFfF = sum(sum(eF2-eF))
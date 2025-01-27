% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%%
% Load the variables from the file into the workspace
%load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\halfwaveDipole.mat');
% f=f0List

dip.pos = [0, 0, 0.1;
           0.1, 0.2, 0.3];
dip.dir = [0, 0, 1;
           1, 0, 0];
dip.complAmpl = [1;2];
f=1e6;

%% quadrature

Nleb = 5810;

[points, weigths, ~] = utilities.getLebedevSphere(Nleb);

%inicialization
construct   = utilities.constants.giveConstants;

%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;
rFar = 1e20/k;
rObserved = points*rFar;
%rObserved = [1, 2, 3]*rFar;
tic
[eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f);
toc


[mF] = fieldEvaluation.magFieldM2(rObserved,dip,f);


tmp = 0.5*real(sum((cross( eF, conj(mF),2).*points ),2));

tic
[fF] = fieldEvaluation.farField(rObserved,dip,f);
propFactor = fieldEvaluation.computePropagationFactor(fF, rObserved, k);
eF2 = fF .* propFactor;
toc

%numerical
integral    = sum( tmp .* weigths);
integral    = integral*rFar^2

result      = fieldEvaluation.powerRadiated(f,dip)
error       = result-integral

integral2   = sum(sum(fF.*conj(fF),2).* weigths)/(2*construct.Z0)
error       = result-integral2
%sum of all errors
sum(sum(eF2-eF))

% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%%
% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\halfwaveDipole.mat');
f=f0List;

%% quadratures
Nleb = 302;
%rObserved = [1, 2, 3]*rFar;
tic
integral = fieldEvaluation.powerQuadrature(Nleb, dip, f)
toc

tic
integral2 = fieldEvaluation.powerQuadratureFar(Nleb, dip, f)
toc
 
%numerical
 
sum(abs(integral-integral2))

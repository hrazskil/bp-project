% Cleaning
clc; % Clear the command window
clear; % Clear all variables from the workspace
%%
% Load the variables from the file into the workspace
load('C:\Users\kilia\Plocha\gitHub\bp-project\graphical test\halfwaveDipole.mat');

% Define the number of theta and phi points
nTh = 60; 
nPh = 120;

% Create linearly spaced vectors
theta = linspace(0, pi, nTh); 
phi = linspace(0, 2*pi, nPh);


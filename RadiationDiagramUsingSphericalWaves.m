%% test projection into spherical waves
clc; clear; close all;

% Load dipole variables from file
load([pwd,'\graphical test\test_structure_1\halfwaveDipole.mat']);

% Define fixed Lebedev quadrature degree (e.g., degree 302 for now)
degree = 302;

% Define physical constants
construct = utilities.constants.giveConstants();
omega = 2 * pi * f0List;  % Angular frequency
k = omega / construct.c0;    % Wavenumber
a = max(sqrt(sum(dip.pos.^2, 2))); % radius of circumscribing sphere
rObs = a*1.5;  % observation radius
ka = k*a; % electrical size

% Generate Lebedev quadrature points and weights for fixed degree
[points, weights, ~] = utilities.getLebedevSphere(degree);
rObserved = points * rObs;  % Scale points to observation distance

% maximum degree of spherical waves
iota = 2;
lmax = ceil(ka + iota*(ka)^(1/3) + 3);

% evaluate electric near field at Lebedev's points
EsCart = fieldEvaluation.eleFieldM2(rObserved, dip, f0List);

% project to spherical waves
fSW = utilities.projectEsTof(lmax, k, EsCart, rObserved, weights);

%% 3D plot of far field

indexVec = ones(size(fSW,1),1);

theta = linspace(0, pi, 60);
phi = linspace(0, 2*pi, 120);

% evaluate far field
farfield = sphericalVectorWaves.sphericalFarField(f0List, fSW, indexVec, theta, phi);

% test evaluation of readiated power
Prad1 = fieldEvaluation.powerQuadratureFar(degree, dip, f0List); % direct evaluation
Prad2 = farfield.Prad; % evaluation from spherical waves
[Prad1; Prad2];

%%%%%%%%%% This does not work without external packages %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hndl = results.plotFarField('farField', abs(farfield.D), ...
%     'theta', theta, 'phi', phi);
% view(42,27);
% 
% export_fig DirectivityFigure.png -m2.5
% 
% hndl.type = 'farField';
% export.plot2TikZ(hndl,...
%    'Azimuth',42,...
%    'Elevation',27,...
%    'Compression',true,...
%    'TikZ_Opacity',1,...
%    'FileName','DirectivityFigure',...
%    'CoordinatesCross',true,...
%    'Colorbar',false);

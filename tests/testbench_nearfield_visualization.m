%% --- Initialization for BTS data without DipoleArray.mat ---
clc; clear; close all;
load('tests/test_structure_3/BTS.mat')
utilities.visualizations.PowerFluxXZ(dip, f0List, [-0.5, 2.4], [-0.5, 0.5], 200);
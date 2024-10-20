function [constStruct] = giveConstants()
constStruct.c0 = 299792458; % speed of light [m/s]
constStruct.eps = 8.854187817e-12; %permitivity of vacuum F·m^(−1)
constStruct.mi = 1/(constStruct.c0^2*constStruct.eps); %permeability of vacuum H·m^(−1)
constStruct.Z0 = sqrt(constStruct.mi/constStruct.eps);
end
function [errorSize] = objectiveFunction(f,dip,dip2,Nleb)
%inicialization
construct   = utilities.constants.giveConstants;
%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;
%% quadrature
Nleb = 302;

% quadrature points and weights
% To integrate a function F(x,y,z) over a unit sphere, call
% sum(F(points).*weigths)
[points, weigths, ~] = utilities.getLebedevSphere(Nleb);

rFar = 1e6/k;
rObserved = points*rFar;

[eF] = fieldEvaluation.eleFieldM2(rObserved,dip,f0List);

[mF] = fieldEvaluation.magFieldM2(rObserved,dip,f0List);

tmp = 0.5*real(sum((cross( eF, conj(mF),2).*points ),2));

[fF] = fieldEvaluation.farField(rObserved,dip,f0List);

% numerical
integral    = sum( tmp .* weigths);
integral    = integral*rFar^2;

result      = fieldEvaluation.powerRadiated(f0List,dip);
errorSize       = result-integral;

sum(sum(fF.*conj(fF),2).* weigths)/(2*construct.Z0);
end


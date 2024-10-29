clear all
close all
clc

% degree: { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,
%   350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,
%   3470, 3890, 4334, 4802, 5294, 5810 };

Nleb = 302; % degree from the list above

% quadrature points and weights
% To integrate a function F(x,y,z) over a unit sphere, call
% sum(F(points).*weigths)
[points, weigths, ~] = utilities.getLebedevSphere(Nleb);

% integrate sin(theta)^2*sin(phi)^2 over a unit sphere
% r = 1

x = points(:,1);
y = points(:,2);
z = points(:,3);

% x.^2 + y.^2 = sin(theta)^2
% y.^2 ./ (x.^2 + y.^2) = sin(phi)^2
% y.^2 = sin(theta)^2*sin(phi)^2

% numerical
sum((y.^2) .* weigths)

% analytical
4*pi/3
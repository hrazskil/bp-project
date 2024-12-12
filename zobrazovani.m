nTh = 51;
nPh = 101;

theta = linspace(0,pi,nTh);
phi = linspace(0,2*pi,nPh);

[Theta, Phi] = meshgrid(theta,phi);

Theta = reshape(Theta,[nTh*nPh,1]);
Phi = reshape(Phi,[nTh*nPh,1]);

[x,y,z] = utilities.transforms.sph2cartKH(ones(nTh*nPh,1),Theta,Phi);
rObserved = [x,y,z];

tic
[fF]    = fieldEvaluation.farField(rObserved,dip,f0List);
toc

[Fr,Fth,Fph] = utilities.transforms.cart2sphVec(x,y,z, fF(:,1), fF(:,2), fF(:,3));

Theta = reshape(Theta, [nPh,nTh]);
Phi = reshape(Phi, [nPh,nTh]);
Fth = reshape(Fth, [nPh,nTh]);
Fph = reshape(Fph, [nPh,nTh]);


figure
contourf(Theta/pi,Phi/(2*pi),real(Fth))
grid on

figure
contourf(Theta/pi,Phi/(2*pi),imag(Fth))
grid on


figure
contourf(Theta/pi,Phi/(2*pi),real(Fph))
grid on

figure
contourf(Theta/pi,Phi/(2*pi),imag(Fph))
grid on

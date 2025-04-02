function giveContourF(dip,f0List,nTh,nPh)

% Create linearly spaced vectors
theta = linspace(0, pi, nTh); 
phi = linspace(0, 2*pi, nPh);

% Create a meshgrid
[Phi, Theta] = meshgrid(phi,theta);

% Reshape into column vectors
Theta = Theta(:); % Flatten into a single column vector
Phi = Phi(:); % <=> Phi = reshape(Phi, [nTh*nPh, 1]);


% Convert spherical to Cartesian
[x, y, z] = utilities.transforms.sph2cartCoor(ones(nTh*nPh, 1), Theta, Phi);
rObserved = [x, y, z]; % Combine Cartesian coordinates into a matrix

% Evaluate the far field
[fF] = fieldEvaluation.farField(rObserved, dip, f0List);

% Transform Cartesian field components to spherical
[~, Fth, Fph] = utilities.transforms.cart2sphVec(...
                x, y, z, fF(:, 1), fF(:, 2), fF(:, 3) ...
                                                 );
% Reshape Theta, Phi, Fth, and Fph for contour plotting
Theta = reshape(Theta, [nTh, nPh]);
Phi = reshape(Phi, [nTh, nPh]);
Fth = reshape(Fth, [nTh, nPh]);
Fph = reshape(Fph, [nTh, nPh]);

% Visualization of Fth and Fph using contour plots with color bars

% Plot for the real part of Fth
figure
contourf(Theta/pi, Phi/(pi), real(Fth));
xlabel('\theta/ \pi')
ylabel('\phi/ \pi')
title('Re[Fth]');
colorbar; % Add color bar for reference


% Plot for the imaginary part of Fth
figure
contourf(Theta/(pi), Phi/pi, imag(Fth));
xlabel('\theta/ \pi')
ylabel('\phi/ \pi')
title('Im[Fth]');
colorbar; 


% Plot for the real part of Fph
figure
contourf(Theta/pi, Phi/(pi), real(Fph));
xlabel('\theta/ \pi')
ylabel('\phi/ \pi')
title('Re[Fph]');
colorbar; 


% Plot for the imaginary part of Fph
figure
contourf(Theta/pi, Phi/(pi), imag(Fph));
xlabel('\theta/ \pi')
ylabel('\phi/ \pi')
title('Im[Fph]');
colorbar;
end
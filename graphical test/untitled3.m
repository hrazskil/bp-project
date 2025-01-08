% Sample data
rObserved = [0, 0, 1; 1, 0, 1; 0, 1, 1]; % Example observed positions
dip.pos = [0, 0, 0]; % Dipole po
dip.complAmpl = 1; % Dipole amplitude
dip.dir = [0, 0, 1]; % Dipole direction
f = 1e9; % Frequency in Hz

% Calculate electric field
eF = fieldevaluation.eleFieldM2(rObserved, dip, f);

% Plot electric field
figure;
quiver3(rObserved(:,1), rObserved(:,2), rObserved(:,3), real(eF(:,1)), real(eF(:,2)), real(eF(:,3)), 'r');
hold on;
quiver3(rObserved(:,1), rObserved(:,2), rObserved(:,3), imag(eF(:,1)), imag(eF(:,2)), imag(eF(:,3)), 'b');
title('Electric Field Visualization');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
legend('Real Part', 'Imaginary Part');
grid on;

% Calculate far-field
fF = farField(rObserved, dip, f);

% Plot far-field radiation pattern
figure;
scatter3(fF(:,1), fF(:,2), fF(:,3), 'filled');
title('Far-Field Radiation Pattern');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;
%%
clc
clear
%%A = rand(numRows, numCols) + 1i * rand(numRows, numCols);
profile on;

profile off;
profile viewer;
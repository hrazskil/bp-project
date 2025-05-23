function PowerFluxXAxis(dip, f, xLim, xPoints, filename)
% Compute and export Poynting vector S_x along the x-axis (y=z=0) for LaTeX plotting

% === Step 1: Define observation points ===
x = linspace(xLim(1), xLim(2), xPoints)';
rObserved = [x, zeros(xPoints,1), zeros(xPoints,1)];

% === Step 2: Evaluate fields ===
E = fieldEvaluation.eleFieldM2(rObserved, dip, f);
H = fieldEvaluation.magFieldM2(rObserved, dip, f);

% === Step 3: Compute S_x ===
S = fieldEvaluation.powerPoynting(E, H);
Sx = S(:,1); % Extract x-component

% === Step 4: Optional clipping ===
Sx(Sx < -250 | Sx > 250) = NaN;

% === Step 5: Export data ===
data = [x, Sx];
writematrix(data, filename, 'Delimiter', '\t'); % saves as .tsv

fprintf('Exported to %s\n', filename);
end

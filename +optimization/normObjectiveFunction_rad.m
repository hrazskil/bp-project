function totalError = normObjectiveFunction_rad(dip, inputData)
%COMPUTERADIATIONERROR Computes the normalized error between
%                      reference and simulated far-field power densities
%                      in vertical and horizontal planes.
%
%   Inputs:
%       dip        - Struct with dipole data, including `complAmpl`
%       inputData  - Struct containing reference power densities, weights,
%                    observation points, and frequency.
%
%   Output:
%       totalError - Combined vertical and horizontal normalized error (scalar)

    % --- Simulate far-field responses on both vertical and horizontal observation planes ---
    farFieldVertical = fieldEvaluation.farFieldM2(inputData.vertical.points, dip, inputData.freq);
    farFieldHorizontal = fieldEvaluation.farFieldM2(inputData.horizontal.points, dip, inputData.freq);

    % --- Compute radiated power densities (|E|^2) ---
    powerDensityVertical = sum(abs(farFieldVertical).^2, 2);
    powerDensityHorizontal = sum(abs(farFieldHorizontal).^2, 2);

    % --- Total radiated power (integrated using weights) ---
    totalPowerVertical = sum(powerDensityVertical .* inputData.vertical.weights);
    totalPowerHorizontal = sum(powerDensityHorizontal .* inputData.horizontal.weights);

    % --- Normalize both simulated and reference radiation data ---
    powerDensityVerticalNorm = powerDensityVertical / totalPowerVertical;
    powerDensityHorizontalNorm = powerDensityHorizontal / totalPowerHorizontal;

    refPowerDensityVerticalNorm = inputData.vertical.rad / inputData.vertical.totalPower;
    refPowerDensityHorizontalNorm = inputData.horizontal.rad / inputData.horizontal.totalPower;

    % --- Compute errors ---
    errorVertical = sum(abs(powerDensityVerticalNorm - refPowerDensityVerticalNorm) .* inputData.vertical.weights);
    errorHorizontal = sum(abs(powerDensityHorizontalNorm - refPowerDensityHorizontalNorm) .* inputData.horizontal.weights);

    % --- Final objective error: total of both planes ---
    totalError = errorVertical + errorHorizontal;
end

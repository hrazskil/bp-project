function giveRadiationPatterns(theta, phi, intensity, frequency)
    % giveRadiationPatterns Plots azimuth and elevation radiation patterns in dB
    % theta: Vector of elevation angles [0, pi]
    % phi: Vector of azimuth angles [0, 2*pi]
    % intensity: Matrix of intensity values (size: [nTheta, nPhi])

    % Find the azimuth and elevation angles for the max intensity
    [intensityMax, maxIntensityIndex] = max(intensity(:));
    [maxThetaIndex, maxPhiIndex] = ind2sub(size(intensity), maxIntensityIndex);

    maxTheta = rad2deg(theta(maxThetaIndex)); % Elevation angle for max intensity
    maxPhi = rad2deg(phi(maxPhiIndex)); % Azimuth angle for max intensity

    % Replace zero intensity values with a small value to avoid log10(0)
    intensity(intensity == 0) = eps;  % eps is a very small value (~2.22e-16)

    % Convert intensity to dB scale
    intensity_dB = 10 * log10(intensity / intensityMax);
    
    %% Create one figure
    figure;

    % Horizontal Plane (Azimuth Pattern) - Looking from Above (XY plane)
    subplot(1, 2, 1); % Left subplot
    thetaIndex = round(size(intensity_dB, 1) / 2); % Pick middle index for θ = 90° (π/2)
    polarplot(phi, intensity_dB(thetaIndex, :), 'b', 'LineWidth', 2);
    title('Azimuth Pattern (XY Plane, θ = 90°)');
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise'); % Adjust direction
    rlim([-40 0]); % Set limits in dB (adjustable)
    grid on;

    % Vertical Plane (Elevation Pattern) - Side View (XZ Plane)
    subplot(1, 2, 2); % Right subplot
    phiIndex = round(size(intensity_dB, 2) / 4); % Pick φ = 0° (assuming φ starts at 0)
    polarplot(theta, intensity_dB(:, phiIndex), 'r', 'LineWidth', 2);
    title('Elevation Pattern (XZ Plane, φ = 0°)');
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    rlim([-40 0]); % Set limits in dB (adjustable)
    grid on;

    %% Bonus Info

    % Calculate the max intensity (in linear scale)
    maxIntensityEng = sprintf('%.2e', intensityMax);

    % Convert frequency to GHz
    f = sprintf('%.2e', frequency / 1e9);

    % Add the bonus info in between the two plots using annotation
    annotation('textbox', [0.3, 0, 0.4, 0.3], 'String', ...
        {...
            sprintf('Frequency [GHz] = %s', f), ...
            sprintf('Max Intensity [W/m2] = %s', maxIntensityEng), ...
            sprintf('Max Gain Azimuth [°]: %.2f', maxPhi), ...
            sprintf('Max Gain Elevation [°]: %.2f', maxTheta) ...
        }, ...
        'edgecolor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

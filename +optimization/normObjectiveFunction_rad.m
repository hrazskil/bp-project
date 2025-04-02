function error = normObjectiveFunction_rad(dip, inputData)
    % Compute the far-field pattern
    fF_vertical = fieldEvaluation.farField(inputData.vertical.points, dip, inputData.freq);
    fF_horizontal = fieldEvaluation.farField(inputData.horizontal.points, dip, inputData.freq);
    
    rad_Vertical = sum(abs(fF_vertical).^2, 2);
    rad_Horizontal = sum(abs(fF_horizontal).^2, 2);
    
    totalPower_ver = trapz(rad_Vertical.*inputData.vertical.weights);
    % sizeStep = inputData.vertical.points(2)-inputData.vertical.points(1); % needs redefining for generalization functions with this examole because 
    % totalPower_ver = sum((rad_Vertical.*inputData.vertical.weights).*sizeStep);
    
    totalPower_hor = trapz(rad_Horizontal.*inputData.vertical.weights);
    % sizeStep = inputData.horizontal.points(2)-inputData.horizontal.points(1); % same as above comment
    % totalPower_hor = sum((rad_Horizontal.*inputData.horizontal.weights).*sizeStep);
    
    % Normalize
    
    rad_Vertical_norm = rad_Vertical / totalPower_ver;
    rad_Horizontal_norm = rad_Horizontal / totalPower_hor;
    
    rad_ref_Vertical_norm = inputData.vertical.rad / inputData.vertical.totalPower;
    rad_ref_Horizontal_norm = inputData.horizontal.rad / inputData.horizontal.totalPower;
    
    % Compute the errors
    error_ver = sum(abs(rad_Vertical_norm - rad_ref_Vertical_norm));       % isn't a virtual zero and is bigger than  error_hor by 1e8
    error_hor = sum(abs(rad_Horizontal_norm - rad_ref_Horizontal_norm));
    
    % Total error as the sum of errors from both planes
    error = error_ver + error_hor;
end
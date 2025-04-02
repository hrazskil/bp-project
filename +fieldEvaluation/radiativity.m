function [radiativity] = radiativity(fF)
construct = utilities.constants.giveConstants();
radiativity = sum(abs(fF).^2, 2)./2*construct.Z0;
end


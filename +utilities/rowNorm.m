function [out] = rowNorm (x)

out = sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);


end
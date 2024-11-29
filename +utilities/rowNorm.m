function [out] = rowNorm (x)
out = sqrt(abs(x(:,1)).^2 + abs(x(:,2)).^2 + abs(x(:,3)).^2);
end
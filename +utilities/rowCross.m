function [C] = rowCross(A, B)
% ROWCROSS Computes the cross product of corresponding rows of two matrices.
%   C = ROWCROSS(A, B) takes two matrices A and B, where each matrix 
%   has three columns representing 3D vectors. The function returns a 
%   matrix C, where each row is the cross product of the corresponding 
%   rows of A and B.
%
%   Inputs:
%       A - An Nx3 matrix where each row is a 3D vector.
%       B - An Nx3 matrix where each row is a 3D vector.
%
%   Outputs:
%       C - An Nx3 matrix containing the cross products of the rows of A and B.


C = [A(:,2).*B(:,3) - A(:,3).*B(:,2), ...
     A(:,3).*B(:,1) - A(:,1).*B(:,3), ...
     A(:,1).*B(:,2) - A(:,2).*B(:,1)];
end


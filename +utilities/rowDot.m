function [D] = rowDot(A, B)
% ROWDOT Computes the dot product of corresponding rows of two matrices.
%   D = ROWDOT(A, B) takes two matrices A and B, where each matrix 
%   has three columns representing 3D vectors. The function returns a 
%   column vector D, where each element is the dot product of the 
%   corresponding rows of A and B.
%
%   Inputs:
%       A - An Nx3 matrix where each row is a 3D vector.
%       B - An Nx3 matrix where each row is a 3D vector.
%
%   Outputs:
%       D - An Nx1 column vector containing the dot products of the 
%           rows of A and B.
D = A(:, 1) .* B(:, 1) + A(:, 2) .* B(:, 2) + A(:, 3) .* B(:, 3);
end



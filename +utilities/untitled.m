% Example vectors
A = [1, 2, 3; 4, 5, 6]; % 2x3 matrix
B = [7, 8, 9; 10, 11, 12]; % 2x3 matrix

% Calculate cross product
C = [A(:,2).*B(:,3) - A(:,3).*B(:,2), ...
     A(:,3).*B(:,1) - A(:,1).*B(:,3), ...
     A(:,1).*B(:,2) - A(:,2).*B(:,1)]
C2 = cross(A,B,2)
% Calculate dot product
D = A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3)

D2 = dot(A,B,2)
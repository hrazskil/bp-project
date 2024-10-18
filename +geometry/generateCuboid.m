function [X,Y,Z]=generateCuboid(depth,length,height,numDepthNodes,numLengthNodes,numHeightNodes,center)

%{
 Inputs:
depth               % length along the x-axis
length
height              % Height of the cylinder
numDepthNodes
numLengthNodes      
numHeightNodes      % Number of nodes along the height
center              % Wanted center of Cuboid coordinates
    % vector giving the Cuboid an orientation
    % Pushing the center of Cuboid to wanted coordinates
    % rotating the Cuboid so it's z-axis goes pararell to a given vector

Output:

%}

arguments
    depth           (1,1) double
    length          (1,1) double
    height          (1,1) double
    numDepthNodes   (1,1) single
    numLengthNodes  (1,1) single
    numHeightNodes  (1,1) single 
    center          (1,3) double
end
% Generate nodal coordinates
x=center(1);
y=center(2);
z=center(3);
x = linspace(-depth+x, depth+x, numDepthNodes);
y = linspace(-length+y, length+y, numLengthNodes);
z = linspace(-height+z, height+z, numHeightNodes);
[X,Y,Z] = meshgrid(x,y,z);
% Reshape the coordinates into column vectors
X = X(:);
Y = Y(:);
Z = Z(:);

%Plot the nodal coordinates
scatter3(X, Y, Z, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Nodal Model of a Cuboid');
axis equal;
end 
function [X,Y,Z]=generateCylinder(rad,height,numCircumNodes,numHeightNodes,center)

%{
    Inputs:
rad             % Radius of the cylinder
height          % Height of the cylinder
numCircumNodes  % Number of nodes along the circumference
numHeightNodes  % Number of nodes along the height
center          % Wanted center of cylinder coordinates
    % vector giving the cylinder an orientation
    % Pushing the center of cylinder to wanted coordinates
    % rotating the cylinder so it's z-axis goes pararell to a given vector

    Output:

%}

arguments
    rad             (1,1) double
    height          (1,1) double
    numCircumNodes  (1,1) single 
    numHeightNodes  (1,1) single 
    center          (1,3) double
end

    % Generate nodal coordinates
x=center(1);
y=center(2);
z=center(3);
phi = 0:2*pi/numCircumNodes:2*pi-2*pi/numCircumNodes;
z = linspace(-height+z, height+z, numHeightNodes);
[phi, Z] = meshgrid(phi, z);
X = rad * cos(phi) + x;
Y = rad * sin(phi) + y;

    % Reshape the coordinates into column vectors
X = X(:);
Y = Y(:);
Z = Z(:);

    %Plot the nodal coordinates
scatter3(X, Y, Z, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Nodal Model of a Cylinder');
axis equal;


end 
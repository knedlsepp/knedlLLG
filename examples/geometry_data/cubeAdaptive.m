function [elements,coordinates] = cubeAdaptive(X,Y,Z,maxIt)
%CUBEADAPTIVE    Generates the mesh of a cube.
%   CUBEADAPTIVE(X,Y,Z,NUMIT) returns a mesh for the cube meshgrid(X,Y,Z)
%   bisecting the simplices NUMIT times adaptively at the minimum
%   Y-coordinate.
%
%   Author: Josef Kemetmueller - 16.12.2013
[X,Y,Z] = meshgrid(X,Y,Z);
coordinates = [X(:), Y(:), Z(:)];
elements = delaunay(coordinates);

[elements,elemGen,coordinates] = validStartMesh(elements,coordinates);

for i = 1:maxIt
    markedCoord = abs((coordinates(:,2)-Y(1)))<eps;
    marked = any(markedCoord(elements),2);
    
    [elements,elemGen,coordinates] = refineC(elements, ...
                                            elemGen,coordinates,marked);
end


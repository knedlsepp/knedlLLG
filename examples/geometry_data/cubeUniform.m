function [elements,coordinates] = cubeUniform(X,Y,Z,numIt)
%CUBEUNIFORM    Generates the mesh of a cube.
%   CUBEUNIFORM(X,Y,Z,NUMIT) returns a mesh for the cube meshgrid(X,Y,Z)
%   bisecting the simplices NUMIT times.
%
%   Author: Josef Kemetmueller - 16.12.2013
percentage = 100/100;
[X,Y,Z] = meshgrid(X,Y,Z);
coordinates = [X(:), Y(:), Z(:)];
elements = delaunay(coordinates);

[elements,elemGen,coordinates] = validStartMesh(elements,coordinates);

for i = 1:numIt
    nE = size(elements,1);
    volumes = getElementVolumes(struct('elements',elements,'coordinates',coordinates));
    [~, I_sort] = sort(volumes);
    marked = I_sort(fix(nE*percentage):end);
    [elements,elemGen,coordinates] = refineC(elements,elemGen,coordinates,marked);
end


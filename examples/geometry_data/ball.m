function [elements,coordinates] = ball(numIt)
%BALL    Generates the approximate mesh of a ball.
%   BALL(NUMIT) returns an approximate mesh of the unit ball bisecting
%   the initial mesh numIt times.
%
%   Author: Josef Kemetmueller - 16.12.2013

%%
percentage = 100/100; % How much to bisect each time.
foo = 2;
switch foo
    case 1
    [X,Y,Z] = meshgrid(-1:2:1,-1:2:1,-1:2:1);
    coordinates = [X(:), Y(:), Z(:)];
    case 2
    coordinates = ...
    [ 0, 0, 0;
     -1, 0, 0;
      1, 0, 0;
      0, 1, 0;
      0,-1, 0;
      0, 0, 1;
      0, 0,-1];
    case 3
    coordinates = ...
    [0.000, 0.000, 1.000 
     0.943, 0.000,-0.333 
    -0.471, 0.816,-0.333 
    -0.471,-0.816,-0.333 ];
end
elements = delaunay(coordinates);

[elements,elemGen,coordinates] = validStartMesh(elements,coordinates);
lengths = @(X) sqrt(sum(X.^2,2));
normalize = @(X) bsxfun(@rdivide, X, lengths(X));


for i = 1:numIt
    nE = size(elements,1);
    volumes = getElementVolumes(struct('elements',elements, ...
                                       'coordinates',coordinates));
    [~, I_sort] = sort(volumes);
    marked = I_sort(fix(nE*percentage):end);
    
    [elements,elemGen,coordinates] = refineC(elements, ...
                                                elemGen,coordinates,marked);
    boundary = getBoundary(elements);
    bNodes = unique(boundary(:));
    coordinates(bNodes,:) = normalize(coordinates(bNodes,:));
end
end

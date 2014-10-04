function plotMagnetization(mesh,mhj,j)
%PLOTMAGNETIZATION    Plots the magnetization.
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   PLOTPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013
scalePlot = abs(max(mesh.coordinates(:,1))-min(mesh.coordinates(:,1)));
scalePlot = 1e-9;
whatToPlot = 'everything';
switch whatToPlot
    case 'boundary'
        frontBd = any(mesh.coordinates>(0.9)*scalePlot/2,2);
        X = @(i) 1/scalePlot*mesh.coordinates(frontBd,i);
        mjhplot = @(i) mhj(frontBd,i);
        quiver3(X(1),X(2),X(3), ...
                mjhplot(1),mjhplot(2),mjhplot(3));
    case 'everything'
        X = @(i) 1/scalePlot*mesh.coordinates(:,i);
        quiver3(X(1),X(2),X(3), ...
                mhj(:,1),mhj(:,2),mhj(:,3));
        xlabel('x');
        ylabel('y');
        zlabel('z');
end

title(sprintf('%d. step',j));
axis vis3d;
axis equal;
drawnow;

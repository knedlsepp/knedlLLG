function plotAndSaveMagnetization(mesh, path, name, mhj, j, variant)
%PLOTANDSAVEMAGNETIZATION    Plots and saves the magnetization.
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   PLOTPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013

plotMagnetization(mesh,mhj,j);
saveMagnetization(mesh,path,name,mhj,j,variant);

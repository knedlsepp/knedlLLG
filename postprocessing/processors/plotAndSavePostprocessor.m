function postProcessor = plotAndSavePostprocessor(mesh, path, name)
%PLOTPOSTPROCESSOR    Constructs a postprocessor that plots and saves.
%   POSTPROCESSOR = PLOTPOSTPROCESSOR(MESH, PATH, NAME) plots the
%   magnetization at the given mesh and saves in PATH under NAME
%
%   See also:
%   PLOTPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013
postProcessor = @(mhj,j,tj) plotAndSaveMagnetization(mesh, path, name, mhj, j, 'mat');
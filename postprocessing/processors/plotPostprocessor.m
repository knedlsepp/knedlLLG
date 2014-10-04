function postProcessor = plotPostprocessor(mesh)
%PLOTPOSTPROCESSOR    Constructs a postprocessor that plots.
%   POSTPROCESSOR = PLOTPOSTPROCESSOR(mesh) plots the magnetization at the
%   given mesh.
%
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013
postProcessor = @(mhj,j,tj) plotMagnetization(mesh,mhj,j);

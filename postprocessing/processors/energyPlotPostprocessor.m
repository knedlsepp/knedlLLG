function postProcessor = energyPlotPostprocessor(mesh, f)
%ENERGYPLOTPOSTPROCESSOR    Constructs a postprocessor that plots the energies.
%   POSTPROCESSOR = ENERGYPLOTPOSTPROCESSOR(MESH, F) plots the energies
%   using the function handle F(t,x) as external field.
%
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   PLOTPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013

postProcessor = @(mhj,j,tj) plotEnergies(mesh, mhj, j, @(x) f(tj,x));
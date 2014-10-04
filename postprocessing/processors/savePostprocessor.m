function postProcessor = savePostprocessor(mesh, path, name)
%SAVEPOSTPROCESSOR    Constructs a postprocessor that saves
%   POSTPROCESSOR = SAVEPOSTPROCESSOR(MESH, PATH, NAME) saves the magnetization
%   in the path PATH using the filename NAME.
%
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   PLOTPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013
postProcessor = @(mhj,j,tj) saveMagnetization(mesh,path,name,mhj,j,'mat');
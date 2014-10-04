kantenlaenge = 39; % [nm]
name = ['standard_problem_3_inhomogen', num2str(kantenlaenge),'nm'];
variant = 'recompute'; % First time always use recompute to build K operator and stuff
%%
thisDir = fileparts(mfilename('fullpath'));
resultFolder = [thisDir,'/../results'];
precomputeFile = [thisDir,'/',name,'_mesh.mat'];
cd(thisDir);
%%
geometry = 'cube24000';
coordinates = dlmread(['coordinates_', geometry, '.txt']);
elements = dlmread(['elements_', geometry, '.txt']);
coordinates = coordinates*kantenlaenge/2*1e-9; % Kantenlaenge 42[nm].
%%
material = struct(...
    'A',  1e-11,...         % [J/m]
    'Js', 1,...              % [T]
    'K',  39788.7346,...     % [J/m^3]
    'easyaxis', [0,0,1]);

%%
switch lower(variant)
    case 'recompute'
        mesh = precomputeMeshValues(coordinates,elements,material);
        save(precomputeFile, 'mesh');
    case 'load'
        load(precomputeFile, 'mesh');
end
%%
t0 = 0;
tEnd = 1e-9;% [ns]
[C_ex, C_ani, tau0, tauEnd, mesh.material.Ms] = nondimensionalization(material, t0, tEnd);
%%
nC = size(coordinates,1);
mh0 = [zeros(nC,1),zeros(nC,1),ones(nC,1)];
% INHOMOGEN:
mh0(coordinates(:,1) >= 0, :) = -mh0(coordinates(:,1) >= 0, :);

f = @(X) zeros(size(X));
ExampleData = struct(...
    'mesh', mesh, ...
    'alpha', 1, ...
    'C_ex', C_ex, ...
    'Piht', @(m) hTildeP0(mesh,strayfield(mesh, m)) ...
                -C_ani*hTildeP1(mesh,anisotropy(m,material)), ...
    'fht', constantInTimeExternalFieldTilde(mesh, f), ...
    'mh0', mh0, ...
    'tau0', tau0, ...
    'tauEnd', tauEnd, ...
    'steps', 2500);
%%
timeStepping(ExampleData, @midpoint_newton, ...
           plotAndSavePostprocessor(mesh, resultFolder, name));        
% TODO: Dimensionen wieder hinzufuegen! M = Ms*m


         
         
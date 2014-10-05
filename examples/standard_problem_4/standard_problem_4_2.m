name = 'standard_problem_4_2';
variant = 'recompute';
%% PERMALLOY
material = struct(...
    'A',  1.3e-11,...	% J/m
    'K',  0.0,...       % J/m^3
    'Ms', 8e5,...        % A/m
    'easyaxis', [1,0,0]);

%%
thisDir = fileparts(mfilename('fullpath'));
precomputeFile = [thisDir,'/',name,'_mesh.mat'];
cd(thisDir);
%%
geometry = 'permalloy29483';
coordinates = dlmread(['coordinates_', geometry, '.txt']);
elements = dlmread(['elements_', geometry, '.txt']);
coordinates = coordinates*1e-9;
try
    load([thisDir,'/../results/m_standard_problem_4_s_state_final.mat'],'mh0');
catch
    disp('Please precompute the initial s-state first.');
end
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
tEnd = 1e-9;
[C_ex, C_ani, tau0, tauEnd, mesh.material.Ms, mu0] = nondimensionalization(material, t0, tEnd);

nC = size(coordinates,1);

fTmu0_inTesla = [-35.5,-6.3,0]*1e-3; %[T]
f_nonDim = fTmu0_inTesla/(material.Ms*mu0);

f_Space = @(X) repmat(f_nonDim, size(X,1) ,1);
ExampleData = struct(...
    'mesh', mesh, ...
    'alpha', 0.02, ...
    'C_ex', C_ex, ...
    'Piht', @(m) hTildeP0(mesh,strayfield(mesh, m)) ...
                -C_ani*hTildeP1(mesh,anisotropy(m,material)), ...
    'fht', constantInTimeExternalFieldTilde(mesh, f_Space), ...
    'mh0', mh0, ...
    'tau0', tau0, ...
    'tauEnd', tauEnd, ...
    'steps', 5000);
%%

timeStepping(ExampleData, @midpoint_newton, ...
             plotAndSavePostprocessor(mesh, [thisDir,'/../results'], name))
         
         
% TODO: Dimensionen wieder hinzufuegen! M = Ms*m


         
         
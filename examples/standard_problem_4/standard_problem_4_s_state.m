name = 'standard_problem_4_s_state';
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
%%
switch lower(variant)
    case 'recompute'
        mesh = precomputeMeshValues(coordinates,elements,material);
        save(precomputeFile, 'mesh');
    case 'load'
        load(precomputeFile, 'mesh');
end
%%
nC = size(coordinates,1);

lengths = @(X) sqrt(sum(X.^2,2));
normalize = @(X) bsxfun(@rdivide, X, lengths(X));
mh0 = normalize([ones(nC,1), zeros(nC,2)]);
t0 = 0;
tEnd = 1e-9;
[C_ex, C_ani, tau0, tauEnd, mesh.material.Ms, mu0] = nondimensionalization(material, t0, tEnd);

f_TimeSpace = @(t,X) 1*(1-(t-t0)/(tauEnd-t0))*ones(size(X));
ExampleData = struct(...
    'mesh', mesh, ...
    'alpha', 0.02, ...
    'C_ex', C_ex, ...
    'Piht', @(m) hTildeP0(mesh,strayfield(mesh, m)) ...
                -C_ani*hTildeP1(mesh,anisotropy(m,material)), ...
    'fht', variableInTimeExternalFieldTilde(mesh, f_TimeSpace, 1), ...
    'mh0', mh0, ...
    'tau0', tau0, ...
    'tauEnd', tauEnd, ...
    'steps', 5000);
%%
timeStepping(ExampleData, @midpoint_newton, ...
             plotAndSavePostprocessor(mesh, [thisDir,'/../results'], name));
%% Now simulate for another nanosecond with external field turned off.
load([thisDir,'/../results/m_standard_problem_4_s_state_ext_',num2str(ExampleData.steps)],'mhj');
ExampleData.mh0 = mhj;
ExampleData.fht = (@(t) 0);
ExampleData.steps = 1000;
save([thisDir,'/../results/m_standard_problem_4_s_state_final'],'mh0');
% TODO: Dimensionen wieder hinzufuegen! M = Ms*m

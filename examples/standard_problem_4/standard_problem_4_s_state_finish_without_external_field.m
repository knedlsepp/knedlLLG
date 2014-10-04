name = 'standard_problem_4_s_state_finalize';
%% PERMALLOY
material = struct(...
    'A',  1.3e-11,...	% J/m
    'K',  0.0,...       % J/m^3
    'Ms', 8e5,...        % A/m
    'easyaxis', [1,0,0]);

%%
thisDir = fileparts(mfilename('fullpath'));
precomputeFile = [thisDir,'/standard_problem_4_s_state_mesh.mat'];
cd(thisDir);
%%
load(precomputeFile, 'mesh');
load([thisDir,'/../results/m_standard_problem_4_s_state_5000'],'mhj');
mh0 = mhj;
%%

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
    'fht', @(t) 0, ...
    'mh0', mh0, ...
    'tau0', tau0, ...
    'tauEnd', tauEnd, ...
    'steps', 1000);
%%

% timeStepping(ExampleData, @midpoint_newton, ...
%              plotAndSavePostprocessor(mesh, [thisDir,'/../results'], name))
timeStepping(ExampleData, @midpoint_newton, ...
             plotAndSavePostprocessor(mesh, [thisDir,'/../results'], name));
load([thisDir,'/../results/m_standard_problem_4_s_state_finalize_2000.mat'],'mhj');
mh0 = mhj;
save([thisDir,'/../results/m_standard_problem_4_s_state_final.mat'],'mh0');
% TODO: Dimensionen wieder hinzufuegen! M = Ms*m


         
         
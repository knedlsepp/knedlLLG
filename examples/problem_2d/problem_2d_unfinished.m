name = 'problem_2d';
variant = 'recompute';
%%
material = struct(...
    'A',  1.3e-11,...          % [J/m]
    'Ms', 7.9577e5,...       % [A/m]
    'K',  39788.7346,...     % [J/m^3]
    'easyaxis', [0,0,1]);

%%
thisDir = fileparts(mfilename('fullpath'));
precomputeFile = [thisDir,'/',name,'_mesh.mat'];
cd(thisDir);
%%
% X = linspace(-20e-9, 20e-9, 2);
% Y = linspace(-20e-9, 20e-9, 2);
% Z = linspace(-20e-9, 20e-9, 2);
% [elements,coordinates] = cubeUniform(X,Y,Z,0);
geometry = 'cube3000';
coordinates = dlmread(['coordinates_', geometry, '.txt']);
elements = dlmread(['elementsSTRAY_', geometry, '.txt']);
coordinates = coordinates*20*1e-9;

X = linspace(-20e-9, 20e-9, 25);
Y = linspace(-20e-9, 20e-9, 25);
[X,Y] = meshgrid(X,Y);
coordinates = [X(:),Y(:)];


coordinates = rand(30,2);

coordinates = coordinates*1e-7;
elements = delaunay(coordinates(:,1),coordinates(:,2));
1

%coordinates = rand(30,3); elements = delaunay(coordinates);
%[elements,coordinates] = cubeUniform([-1e9,1e9],[-1e9,1e9],[-1e9,1e9],4);


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
tEnd = 2e-9;
[C_ex, C_ani, tau0, tauEnd] = nondimensionalization(material, t0, tEnd);

nC = size(coordinates,1);
mh0 = [zeros(nC,2),ones(nC,1)];
mh0 = [zeros(nC,1),-ones(nC,1),zeros(nC,1)];



%mh0 = normalize(ones(nC,3));

mh0(coordinates(:,1) >= 0, :) = -mh0(coordinates(:,1) >= 0, [1,3,2]);


lengths = @(X) sqrt(sum(X.^2,2));
normalize = @(X) bsxfun(@rdivide, X, lengths(X));



%mh0 = normalize(ones(nC,3));

mh0 = normalize(1/2-rand(size(coordinates,1),3));

% P = mesh.coordinates;
% fn = mh0;
% quiver3(P(:,1),P(:,2),P(:,3),fn(:,1),fn(:,2),fn(:,3),0.5, 'color','b');


%f_Space = @(X) 1e0*[zeros(size(X,1),1),zeros(size(X,1),1),zeros(size(X,1),1)];
f_Space = @(X) 1e0*[ones(size(X,1),1),zeros(size(X,1),1),zeros(size(X,1),1)];
f_TimeSpace = @(t,X) (tauEnd-t)/(tauEnd-tau0)*f_Space(X);

ExampleData = struct(...
    'mesh', mesh, ...
    'alpha', 0, ...
    'C_ex', C_ex, ...
    ...%'Piht', @(m) hTildeP0(mesh,strayfield(mesh, m)) ...
    ...%            -C_ani*hTildeP1(mesh,anisotropy(m,material)), ...
    'Piht', @(m) zeros(size(m)), ...
    'fht', constantInTimeExternalFieldTilde(mesh, f_Space), ...
    ...%'fht', @(t) hTildeL2(mesh, @(x) f_TimeSpace(t,x), 3), ...
    'mh0', mh0, ...
    'tau0', tau0, ...
    'tauEnd', tauEnd, ...
    'steps', 10000);
%%
%timeStepping(ExampleData, @midpoint_newton, ...
%            plotAndSavePostprocessor(mesh, [thisDir,'/../results'], name))
timeStepping(ExampleData, @midpoint_newton, ...
              energyPlotPostprocessor(mesh, f_TimeSpace))
         
         
% TODO: Dimensionen wieder hinzufuegen! M = Ms*m


         
         
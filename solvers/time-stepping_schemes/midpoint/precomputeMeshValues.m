function mesh = precomputeMeshValues(coordinates, elements, material)
%PRECOMPUTEMESHVALUES    Precomputes data which is often needed in the simulation
%   MESH = PRECOMPUTEMESHVALUES(COORDINATES, ELEMENTS, MATERIAL) returns
%   a mesh structure containing often used data of the mesh.
%   COORDINATES and ELEMENTS represent the mesh by standard simplex-vertex
%   format. MATERIAL is a struct containing the material parameters     
%       'A'  [J/m]      used for exchange energy computation
%       'Ms' [A/m]      used for external and strayfield energy computation
%       'K'  [J/m^3]    used for anisotropy energy computation
%       'easyaxis' []   used for anisotropy energy computation and anistropy.
%
%   Author: Josef Kemetmueller - 16.12.2013
%% Meshdata
mesh.coordinates = coordinates;
mesh.elements = elements;
mesh.volumes = getElementVolumes(mesh);
mesh.betas = getBetas(mesh);
mesh.hatGrads = getHatGrads(mesh);
%%
lengths = @(X) sum(X.^2);
if exist('material','var')
    material.easyaxis = reshape(material.easyaxis,3,[]);
    assert(abs(lengths(material.easyaxis)-1)<1e2*eps, 'Easyaxis not normed');
    mesh.material = material;
    [~, ~, ~, ~, mesh.material.Ms] = nondimensionalization(material, 0, 0);
end
%% inner_gradHatI_gradHatJ
% Used in: inner_HatIxlaphHatJ_P1_h
mesh.inner_gradHatI_gradHatJ = cell(dimMesh(mesh)+1,dimMesh(mesh)+1);
for i = 1:dimMesh(mesh)+1
    for j = 1:dimMesh(mesh)+1
        mesh.inner_gradHatI_gradHatJ{i,j} = ...
            inner_P0_P0_L2(mesh, ...
                           mesh.hatGrads{i},mesh.hatGrads{j},...
                           'elementwise');
    end
end
%% Stiffness matrix
mesh.stiffness = inner_gradP1_gradP1_L2(mesh);
%% Boundary-data and Double-Layer-Potential
% outwards orientation of the boundary-faces is important for buildK.
mesh.bd = struct('elements', getBoundary(mesh,'outwards'), ...
                 'coordinates', coordinates);
% Remove unused nodes from boundary-mesh and renumber.
[mesh.bd, ~, mesh.bd.nodesOf3DMesh] = cleanMesh(mesh.bd);
mesh.bd.volumes = getElementVolumes(mesh.bd); % volumes = areas
% volumes are actually the face areas, but we call them volumes to be
% consistent with the mesh-data-structure.
if dimMesh(mesh)==3
    assert(exist('buildK','file')==3, 'You have to compile the K-operator first. Goto folder external_codes compileMayrKOperator.'); 
    disp('Computing the double layer potential. You can go for a coffee break now');
    mesh.bd.DLPot = bsxfun(@rdivide, -buildK(mesh.bd.coordinates, ...
                                        mesh.bd.elements, mesh.bd.volumes, 0.5),...
                           mesh.bd.volumes);
% The above computation yields the same DLPot independent of the scaling of
% mesh.bd.coordinates. (If you divide mesh.bd.coordinates by 10 and the
% areas by 100 it will return the same matrix)
else
    % Also outwards oriented boundary doesn't work with getBoundary.
    mesh.coordinates = [coordinates, zeros(numCoordinates(mesh), 3-dimSpace(mesh))];
    warning('No strayfield computation is available for %dD', dimMesh(mesh));
    mesh.bd.DLPot = zeros([size(mesh.bd.elements,1),size(mesh.bd.coordinates,1)]);
end
end



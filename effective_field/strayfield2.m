function minusGRADu = strayfield2(mesh, m)
%STRAYFIELD    Computes the strayfield.
%   MINUSGRADU = STRAYFIELD(MESH, M) returns the strayfield computed via
%   the Fredkin-Koehler FEM-BEM-coupling ansatz. The result is a P0
%   function, namely minus the gradient of the magnetostatic potential U.
%   [-gradP1(MESH, U)]
%   
%   Author: Josef Kemetmueller - 16.12.2013
nC = size(mesh.coordinates,1);
dimMesh = size(mesh.elements,2)-1; % Should always be 3...
h = mesh.volumes(1)^(1/dimMesh);
%% STEP 1: compute (\nabla u1, \nabla v_h) + h^(-2)*(u1,1)(v_h,1) = (M,\nabla v_h)
%  ------  for all v_h in S^1_*(T_h)
A = mesh.stiffness + h^(-4)*mesh.betas*mesh.betas';
% assembe right-hand side via (grad(V_j),m)_L2
b = inner_GradHatI_P1_L2(mesh, m);
x1 = A \ b;
%% STEP 2: compute Dirichlet data g_h = J_h(K - 1/2)u1|_Gamma
%  ------
u1 = x1(mesh.bd.nodesOf3DMesh);
% use 'double layer' potential operator K-1/2
gh_P0 = mesh.bd.DLPot*u1;
% get P1 approximation by Clement-interpolation on boundary
g_h = interpolateClementP0(mesh.bd,gh_P0);

%% STEP 3: compute solution of (\nabla u2h, \nabla vh) = 0 
%  ------  with u2h|_Gamma = g_h
% prescribe values at boundary
x2 = zeros(nC,1);
x2(mesh.bd.nodesOf3DMesh) = g_h;
% assemble right-hand side
b = -mesh.stiffness*x2;
% compute P1-FEM solution
freenodes = setdiff(1:nC, mesh.bd.nodesOf3DMesh);
x2(freenodes) = mesh.stiffness(freenodes, freenodes)\b(freenodes);

%% STEP 4: obtain discrete solution and demagnetization field
%  ------
% compute demagnetization field 
% PhM = \grad x1 + \grad x2 on elements

minusGRADu = -gradP1(mesh, x1+x2);
function minusGRADu = strayfield(mesh, m)
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
%% STEP 1: compute (\nabla u1, \nabla v_h) = (M,\nabla v_h)
%  ------  for all v_h in S^1_*(T_h)
% compute P1-FEM approx. of Neumann problem by solving extended system,
% where the extra row and column results in a lagrange multiplier, so for
% the solution holds \int_(T_h) x1 = 0.
% We scale the lagrange multiplier, so the matrix-condition is ok.
% tic;
% EA = [mesh.stiffness, h^(-2)*mesh.betas; h^(-2)*mesh.betas' 0];
% % assembe right-hand side via (grad(V_j),m)_L2
% Eb = [inner_GradHatI_P1_L2(mesh, m); 0];
% Ex1 = EA \ Eb;
% x1 = Ex1(1:end-1,1);
% toc;
%%
internalNodes = setdiff(mesh.elements(:),mesh.bd.nodesOf3DMesh);
free = (1:nC)~=internalNodes(1);
x1 = zeros(nC,1);
b = inner_GradHatI_P1_L2(mesh, m);
x1(free) = mesh.stiffness(free,free)\b(free);
%%
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
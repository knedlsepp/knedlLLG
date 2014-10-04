function laph = laplace_hP1(mesh, P1)
%LAPLACE_hP1 Discrete version of the laplacian.
%   LAPH = LAPLACE_hP1(MESH, P1) returns the discrete (mass lumped)
%   laplacian of the piecewise linear, globally continuous function P1.
%   The result satisfies:
%   (LAPH, P1s)_h = -(grad P1, grad P1s)_L2 for all piecewise linear, globally
%   continuous P1s.
%
%   Author: Josef Kemetmueller - 16.12.2013
laph = -bsxfun(@rdivide, mesh.stiffness*P1, mesh.betas);
% This would be an alternative to compute it, but since its a diagonal
% matrix, we do it explicitly:
%laph = -inner_P1_P1_h(mesh)\inner_gradP1_gradP1_L2(mesh,P1);
end
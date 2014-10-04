function energyEx = energyExchange(mesh, m)
%ENERGYEXCHANGE    Compute the exchange energy
%   ENERGYEXCHANGE(MESH, M) returns the exchange energy of the
%   magnetziation M.
%
%   Author: Josef Kemetmueller - 16.12.2013
energyEx = mesh.material.A*inner_gradP1_gradP1_L2(mesh, m, m);
% Alternatively:
% energyEx = mesh.material.A*inner_P0_P0_L2(mesh, gradP1(mesh, m), gradP1(mesh, m));
function energyStray = energyStrayfield(mesh, m)
%ENERGYSTRAYFIELD    Compute the strayfield energy
%   ENERGYSTRAYFIELD(MESH, M) returns the strayfield energy of the
%   magnetziation M.
%
%   Author: Josef Kemetmueller - 16.12.2013
[~, ~, ~, ~, mesh.material.Ms] = nondimensionalization(mesh.material, 0, 0);


mu0 = 4*pi*1e-7;    % vacuum permeability in [T*m/A]
SF = strayfield(mesh, m);
energyStray = -mu0*mesh.material.Ms^2/2*inner_P0_P1_L2(mesh,SF,m);
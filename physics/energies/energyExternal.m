function energyExt = energyExternal(mesh, m, f_t, quadDeg)
%ENERGYEXTERNAL    Compute the external [Zeeman] energy
%   ENERGYEXTERNAL(MESH, M, F_T, QUADDEG) returns the external energy of
%   the magnetziation M. F_T(X) is a function handle for the external field
%   at some time T and QUADDEG is the quadrature degree.
%
%   Author: Josef Kemetmueller - 16.12.2013

[~, ~, ~, ~, mesh.material.Ms] = nondimensionalization(mesh.material, 0, 0);
mu0 = 4*pi*1e-7; % vacuum permeability in [T*m/A]     
energyExt = -mu0*mesh.material.Ms^2*inner_L2_P1_L2(mesh, f_t, m, quadDeg);
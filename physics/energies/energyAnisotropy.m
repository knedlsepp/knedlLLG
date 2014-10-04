function energyAni = energyAnisotropy(mesh, m)
%ENERGYANISOTROPY    Compute the anisotropy energy
%   ENERGYANISOTROPY(MESH, M) returns the anisotropy energy of the
%   magnetziation M.
%
%   Author: Josef Kemetmueller - 16.12.2013
energyAni = mesh.material.K*inner_P1_P1_L2(mesh, 1-m*mesh.material.easyaxis, ...
                                                 1+m*mesh.material.easyaxis);
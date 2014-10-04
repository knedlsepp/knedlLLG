function M = inner_P1xlaphHatJ_HatI_h(mesh, P1)
%INNER_P1XLAPHHATJ_HATI_h    Laplace-Cross product matrix using mass lumping.
%   M = INNER_P1XLAPHHATJ_HATI_h(MESH, P1) returns the matrix defined by:
%   M(i,j) = (P1 x lap_h(hat(j)) , hat(i))_h, where the hat functions are
%   numbered by hat(nC*(dim-1)+node), dim=1...3, node=1...nC. 
%
%   Author: Josef Kemetmueller - 16.12.2013
M = -inner_HatIxlaphHatJ_P1_h(mesh, P1);
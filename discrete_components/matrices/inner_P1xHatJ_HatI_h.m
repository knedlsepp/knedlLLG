function M = inner_P1xHatJ_HatI_h(mesh, P1)
%INNER_P1XHATJ_HATI_H    Cross product matrix using mass lumping.
%   M = INNER_P1XHATJ_HATI_H(MESH, P1) returns the matrix defined by:
%   M(i,j) = (P1 x hat(j) , hat(i))_h, where the hat functions are numbered
%   by hat(nC*(dim-1)+node), dim=1...3, node=1...nC.
%   This results in the following block matrix, which is skew-symmetric.
%
%       [(m(n)xe1,e1)_h, (m(n)xe2,e1)_h, (m(n)xe3,e1)_h ]
%   M = [(m(n)xe1,e2)_h, (m(n)xe2,e2)_h, (m(n)xe3,e2)_h ]
%       [(m(n)xe1,e3)_h, (m(n)xe2,e3)_h, (m(n)xe3,e3)_h ]
%
%   Author: Josef Kemetmueller - 16.12.2013
M = -inner_HatIxHatJ_P1_h(mesh, P1);
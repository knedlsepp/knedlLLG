function M = inner_HatIxHatJ_P1_h(mesh, P1)
%INNER_HATIXHATJ_P1_H    Cross product matrix using mass lumping.
%   M = INNER_HATIXHATJ_P1_H(MESH, P1) returns the matrix defined by:
%   M(i,j) = (hat(i) x hat(j) , P1)_h, where the hat functions are numbered
%   by hat(nC*(dim-1)+node), dim=1...3, node=1...nC.
%   This results in the following block matrix, which is skew-symmetric.
%
%       [(e1xe1,m(n))_h, (e1xe2,m(n))_h, (e1xe3,m(n))_h ]
%   M = [(e2xe1,m(n))_h, (e2xe2,m(n))_h, (e2xe3,m(n))_h ]
%       [(e3xe1,m(n))_h, (e3xe2,m(n))_h, (e3xe3,m(n))_h ]
%
%   Author: Josef Kemetmueller - 16.12.2013
assert(all(size(P1) == size(mesh.coordinates)),'Input not nC-by-3 function.');
%%
nC = numCoordinates(mesh);
dim = 3;
%%
indLocal2Global = @(d,nodes) bsxfun(@plus,nC*(d-1),reshape(nodes,[],1));
% Yields the indices to the global block-matrix: <e(i) x e(j) , phi>
ILocal2Global = @(nodes) indLocal2Global([1 1 2 2 3 3],nodes);
JLocal2Global = @(nodes) indLocal2Global([2 3 1 3 1 2],nodes);
% Yields the entries of the block-matrix: <e(i) x e(j) , phi>
crossMatLocal2Global = @(phi) [           phi(:,3) -phi(:,2) ...
                               -phi(:,3)            phi(:,1) ...
                                phi(:,2) -phi(:,1)          ];
% The non-mass-lumped version would probably be a bit more complicated.
S = crossMatLocal2Global(inner_P1_P1_h(mesh,P1));
I = ILocal2Global(1:nC); 
J = JLocal2Global(1:nC);

M = sparse(I,J,S,dim*nC,dim*nC);
end
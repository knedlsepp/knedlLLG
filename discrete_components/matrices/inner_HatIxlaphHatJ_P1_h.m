function M = inner_HatIxlaphHatJ_P1_h(mesh, P1)
%INNER_HATIXLAPHHATJ_P1_h    Laplace-Cross product matrix using mass lumping.
%   M = INNER_HATIXLAPHHATJ_P1_h(MESH, P1) returns the matrix defined by:
%   M(i,j) = (hat(i) x lap_h(hat(j)) , P1)_h, where the hat functions are
%   numbered by hat(nC*(dim-1)+node), dim=1...3, node=1...nC. 
%
%   Author: Josef Kemetmueller - 16.12.2013
assert(all(size(P1) == size(mesh.coordinates)), 'Input not nC-by-3 function.');
%%
nC = numCoordinates(mesh);   % number of nodes
dim = 3;
indLocal2Global = @(d,nodes) bsxfun(@plus,nC*(d-1),reshape(nodes,[],1));
% Yields the indices to the global block-matrix: <e(i) x e(j) , P1>
ILocal2Global = @(nodes) indLocal2Global([1 1 2 2 3 3],nodes);
JLocal2Global = @(nodes) indLocal2Global([2 3 1 3 1 2],nodes);
% Yields the entries of the block-matrix: <e(i) x e(j) , P1>
crossMatLocal2Global = @(X) [         X(:,3) -X(:,2) ...
                             -X(:,3)          X(:,1) ...
                              X(:,2) -X(:,1)        ];
[S,I,J] = deal(zeros(numElements(mesh),6,dimMesh(mesh)+1,dimMesh(mesh)+1));

for i = 1:dimMesh(mesh)+1
    phiElI = P1(mesh.elements(:,i),:);
    for j = 1:dimMesh(mesh)+1
        S(:,:,i,j) = bsxfun(@times,-crossMatLocal2Global(phiElI), ...
                                   mesh.inner_gradHatI_gradHatJ{i,j});
        I(:,:,i,j) = ILocal2Global(mesh.elements(:,i));
        J(:,:,i,j) = JLocal2Global(mesh.elements(:,j));
    end
end
M = sparse(I(:),J(:),S(:),dim*nC,dim*nC);
end
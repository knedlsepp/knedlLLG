function P1 = hTildeL2(mesh, L2, quadDeg)
%HTILDEL2    'Isometry' from ({P1}, (.,{P1})_L2) to ({P1}, (.,{P1})_h)
%   P1 = hTildeL2(mesh, L2) returns the elementwise linear, globally
%   continuous function P1, that satisfies:
%   $(P1, P1s)_h = (L2, P1s)_L2$ for all elementwise linear, globally
%   continuous P1s.
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('quaddeg','var');
    quadDeg = 3;
end
assert(size(L2(mesh.coordinates(1,:)),1)==1,'L2 values must be row vectors.');

SPs = inner_L2_P1_L2(mesh, L2, [], quadDeg);
P1 = bsxfun(@rdivide, SPs, getBetas(mesh));

% This would be an alternative to compute it, but since its a diagonal
% matrix, we do it explicitly:
%P1 = inner_P1_P1_h(mesh)\inner_L2_P1_L2(mesh, L2, [], quadDeg);
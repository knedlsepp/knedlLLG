function P1_h = hTildeP1(mesh, P1_L2)
%HTILDEP1    'Isometry' from ({P1}, (.,{P1})_L2) to ({P1}, (.,{P1})_h)
%   P1_h = hTildeP1(mesh, P1_L2) returns the elementwise linear, globally
%   continuous function P1_h, that satisfies:
%   $(P1_h, P1s)_h = (P1_L2, P1s)_L2$ for all elementwise linear, globally
%   continuous P1s.
%
%   Author: Josef Kemetmueller - 16.12.2013

P1_h = bsxfun(@rdivide, inner_P1_P1_L2(mesh, P1_L2), getBetas(mesh));

% This would be an alternative to compute it, but since its a diagonal
% matrix, we do it explicitly:
%P1 = inner_P1_P1_h(mesh)\inner_P1_P1_L2(mesh, P1_L2);

end
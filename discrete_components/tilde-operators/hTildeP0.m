function P1 = hTildeP0(mesh, P0)
%HTILDEP0    'Isometry' from ({P0}, (.,{P1})_L2) to ({P1}, (.,{P1})_h)
%   P1 = hTildeP0(mesh, P0) returns the elementwise linear, globally
%   continuous function P1, that satisfies:
%   $(P1, P1s)_h = (P0, P1s)_L2$ for all elementwise linear, globally
%   continuous P1s.
%
%   Author: Josef Kemetmueller - 16.12.2013

SPs = inner_P0_P1_L2(mesh, P0);
P1 = bsxfun(@rdivide, SPs, getBetas(mesh));

% This would be an alternative to compute it, but since its a diagonal
% matrix, we do it explicitly:
%P1 = inner_P1_P1_h(mesh)\inner_P0_P1_L2(mesh, P0);
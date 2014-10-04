function f = variableInTimeExternalFieldTilde(mesh, f_TimeSpace, quadDeg)
%VARIABLEINTIMEEXTERNALFIELDTILDE    Generate function handle which is
%variable in time.
%   F = VARIABLEINTIMEEXTERNALFIELDTILDE(MESH, F_L2, QUADDEG) returns a
%   function handle f(t), which yields a discrete P1 field for each
%   time t. The input should be a function handle F_L2(t,X) which for every
%   time t can be evaluated using an array whose rows are points in R^3.
%
%   Author: Josef Kemetmueller - 16.12.2013

f = @(t) hTildeL2(mesh, @(X) f_TimeSpace(t,X), quadDeg);
end
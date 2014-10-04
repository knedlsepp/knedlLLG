function f = constantInTimeExternalFieldTilde(mesh, f_Space, quadDeg)
%CONSTANTINTIMEEXTERNALFIELDTILDE    Generate function handle which is
%constant in time.
%   F = CONSTANTINTIMEEXTERNALFIELDTILDE(MESH, F_SPACE, QUADDEG) returns a
%   function handle f(t), which yields a discrete P1 field for each
%   time t. The input should be a function handle F_SPACE(X) which can be
%   evaluated using an array whose rows are points in R^3.
%
%   Author: Josef Kemetmueller - 16.12.2013

if ~exist('quadDeg','var')
    quadDeg = 3;
end
timeConstant = hTildeL2(mesh, f_Space, quadDeg);
f = @(t) timeConstant;
end
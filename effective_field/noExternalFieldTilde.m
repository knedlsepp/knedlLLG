function f = noExternalFieldTilde()
%NOEXTERNALFIELDTILDE    Generate function handle for no external field.
%   F = NOEXTERNALFIELDTILDE() returns a
%   function handle f(t), which yields a zero-P1 field for each
%   time t. 
%
%   Author: Josef Kemetmueller - 16.12.2013

f = @(t) 0;
% Theoretically we would need a zero-Matrix, but this works too.
end
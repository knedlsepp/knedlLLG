function [C_ex, C_ani, tau0, tauEnd, Ms, mu0] = nondimensionalization(material, t0, tEnd)
%NONDIMENSIONALIZATION    Computes the constants necessary for the
%nondimensional formulation of the LLG.
%
%   Author: Josef Kemetmueller - 16.12.2013

% physical constants
gamma_0 = 2.210173e5;    % gyromagnetic ratio in [m/(As)]
mu0 = 4*pi*1e-7;    % permeability of vacuum in [T*m/A]     
L = 1;% [m]
assert(isfield(material,'Js')||isfield(material,'Ms'), 'Either Js or Ms must be given.');

% saturation magnetization Ms in [A/m]
if ~isfield(material,'Ms')
    Ms = material.Js/mu0;
else
    Ms = material.Ms;
end
C_ex  = 2*material.A/(mu0*Ms*Ms*L*L);
C_ani = material.K/(mu0*Ms*Ms);
tau0 = gamma_0*Ms*t0;
tauEnd = gamma_0*Ms*tEnd;
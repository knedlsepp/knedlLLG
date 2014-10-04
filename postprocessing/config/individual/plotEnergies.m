function plotEnergies(mesh, mhj, j, f_t)
%PLOTENERGIES    Plots the energies.
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   PLOTPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013
numSteps = 10;
persistent energies;
if j==1
    energies = struct('Exchange',[],...
                      'Anisotropy', [], ...
                      'Strayfield', [], ...
                      'External',[]);
end

energyExch = energyExchange(mesh, mhj);
energyAni = energyAnisotropy(mesh, mhj);
energyStray = energyStrayfield(mesh, mhj);
energyExt = energyExternal(mesh, mhj, f_t, 2);

energies.Exchange(end+1) = energyExch;
energies.Anisotropy(end+1) = energyAni;
energies.Strayfield(end+1) = energyStray;
energies.External(end+1) = energyExt;
Esum = energies.Exchange+...
       energies.Anisotropy+...
       energies.Strayfield+...
       energies.External;

% if mod(j,numSteps)~=0
%     return;
% end
figure(1);
plotMagnetization(mesh,mhj,j)

figure(2);
clf;
plot([energies.Exchange; energies.Anisotropy; ...
      energies.Strayfield; energies.External; Esum]');
legend('Exchange', 'Anisotropy', 'Strayfield', 'External', 'Sum of energies');
drawnow;

end
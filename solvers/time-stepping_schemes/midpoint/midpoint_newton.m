function mhjpe = midpoint_newton(mesh, alpha, C_ex, Piht, fht, mhj, tj, tjpe)
%MIDPOINT_NEWTON    Computes the next magnetization using the midpoint
%scheme and a newton-iteration to solve the nonlinear system.
%   MHJPE = MIDPOINT_NEWTON(MESH, ALPHA, C_EX, PIHT, FHT, MHJ, TJ, TJPE)
%   returns the magnetization m_{j+1} for the given input
%       mesh         Geometry of the problem
%       alpha        damping factor
%       C_ex         Exchange constant
%       @(X) Piht(X) Discretized Pi-operator
%       @(t) fht(X)  Discretized external field
%       mhj          Magnetization at previous timestep tj
%       tj           Previous time
%       tjpe         Next time    
%
%   Author: Josef Kemetmueller - 16.12.2013

maxIt = 40;
k = tjpe-tj;
%%
h = norm(diff(mesh.coordinates(mesh.elements(1,1:2),:),1));
epsilon = 1e-9*k*(h^2); % TODO: What is the correct epsilon?
%%
mhjpeellm1 = mhj;
%%
hack = false;
if ~hack
h_eff = @(mhj, mhjpe) C_ex*laplace_hP1(mesh,(mhjpe+mhj)/2) ...
                      + Piht(mhj) + fht((tjpe+tj)/2);
else
    warning('implicit-hack.');
h_eff = @(mhj, mhjpe) C_ex*laplace_hP1(mesh,(mhjpe+mhj)/2) ...
                      + Piht((mhj+mhjpe)/2) + fht((tjpe+tj)/2);
end
F = @(mhj, mhjpe) inner_P1_P1_h(mesh, [],...
                  1/k*(mhjpe - mhj) ...
                  - alpha*cross(mhj, 1/k*(mhjpe - mhj),2) ...
                  + cross((mhjpe+mhj)/2, h_eff(mhj,mhjpe),2));
change(1) = NaN;
for l = 1:maxIt
    dF = 1/k*inner_P1_P1_h(mesh,[],[],3) ...
        + inner_P1xHatJ_HatI_h(mesh, -alpha/k*mhj ...
                                    -1/2*(h_eff(mhj,mhjpeellm1))) ...
        + C_ex/2*inner_P1xlaphHatJ_HatI_h(mesh,(mhj+mhjpeellm1)/2);
    % SOLVE
    tmp = reshape(F(mhj,mhjpeellm1),[],1); %HACK: Just to see profiling speed
    mhjpeell = mhjpeellm1 - reshape(dF\tmp,[],3);    
    %% Stopping criterion
    change(l+1) = norm_P1_h(mesh, mhjpeellm1-mhjpeell);
    if  (change(l+1)<epsilon)
        fprintf('Timestep completed after %d steps.\n',l);
        break;
    elseif change(l+1)*(1+1e-8) > change(l)
        warning(['Newton tolerance not reached, but the iteration didn''t', ...
                 'change anymore.\n We would have wished for a tolerance ', ...
                 'of eps=%g, but only reached eps=%g'], epsilon, change(end));
        fprintf('Here is a table of the changes in norm: \n');
        fprintf('Step %d: %g\n', [1:length(change); change]);
        break;
    else % Continue iterating.
        mhjpeellm1 = mhjpeell;
    end
end
assert(l<maxIt,'Maximum number of iterations reached.');
mhjpe = mhjpeell;
end
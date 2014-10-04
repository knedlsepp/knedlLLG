function timeStepping(ExampleData, timestepFunction, postprocessMagnetization)
%TIMESTEPPING    Solves the LLG-equation using uniform timestepping.
%   TIMESTEPPING(EXAMPLEDATA, TIMESTEPFUNCTION, POSTPROCESSOR) performs
%   uniform timestepping to solve the nondimensional LLG-equation:
%  --------------------------------------------------
%       m_t - alpha*cross(m,m_t) = -cross(m,h_eff)
%           where h_eff is defined by:
%       h_eff := C_ex*laplace(m) + Pi(m) + f(t)
%       Pi(m) := -C_ani*DPhi(m) + strayfield(m)
%  --------------------------------------------------
%   EXAMPLEDATA is a struct, that contains the following data:
%       mesh         Geometry of the problem
%       alpha        damping factor
%       C_ex         Exchange constant
%       @(X) Piht(X) Discretized Pi-operator
%       @(t) fht(X)  Discretized external field
%       tau0         Start time
%       tauEnd       End time
%       steps        Number of timesteps to perform
%
%   TIMESTEPFUNCTION is a function handle to a solver that computes the
%   next timestep. At the moment the options are
%       @midpoint_newton      A newton iteration of the midpoint scheme
%       @midpoint_fixedPoint  A fixedpoint iteration of the midpoint scheme
%   
%   POSTPROCESSOR is function handle, that specifies what to do with each
%   simulation. At the moment the options are:
%       @savePostprocessor(path, name)
%       @plotPostprocessor(mesh)
%       @plotAndSavePostprocessor(mesh, path, name)
%       @energyPlotPostprocessor(mesh, f)
%
%   Author: Josef Kemetmueller - 16.12.2013

mhj = ExampleData.mh0;

t = linspace(ExampleData.tau0, ExampleData.tauEnd, ExampleData.steps+1);
for j = 1:length(t)-1
    mhj = timestepFunction(ExampleData.mesh, ExampleData.alpha, ExampleData.C_ex, ...
                           ExampleData.Piht, ExampleData.fht, mhj, ...
                           t(j), t(j+1));
    postprocessMagnetization(mhj,j,t(j));
end
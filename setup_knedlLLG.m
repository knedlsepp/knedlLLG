% THIS SCRIPT COMPILES THE BISECTION-SOURCE
[pathstrmain, name, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(pathstrmain),'-begin');
setup_p1afem_nD;
compileKOperator;
cd(pathstrmain);

% THIS SCRIPT COMPILES THE BISECTION-SOURCE
[knedlpathstrmain, name, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(knedlpathstrmain),'-begin');
setup_p1afem_nD;
compileKOperator;
cd(knedlpathstrmain);

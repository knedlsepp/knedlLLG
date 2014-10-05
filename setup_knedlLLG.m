% THIS SCRIPT COMPILES THE BISECTION-SOURCE
[knedlpathstrmain, name, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(knedlpathstrmain),'-begin');
assert(exist('setup_p1afem_nD', 'file')==2,'Submodule p1afem-nD not found');
assert(exist('compileKOperator', 'file')==2,'Submodule K-operator not found');
setup_p1afem_nD;
compileKOperator;
cd(knedlpathstrmain);

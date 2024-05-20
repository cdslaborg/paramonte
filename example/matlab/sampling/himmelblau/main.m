warning off;
diary main.out.m;
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

sampler = pm.sampling.Paradram();
sampler.spec.randomSeed = 284651; % make sampling reproducible.
sampler.spec.outputChainSize = 10000; % Use a small chain size for illustration.
sampler.spec.parallelismNumThread = []; % Use these many parallel threads.
sampler.run(@getLogFunc, 2);

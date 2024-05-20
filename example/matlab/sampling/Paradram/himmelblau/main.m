warning off;
diary main.out.m;
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

sampler = pm.sampling.Paradram();
sampler.spec.randomSeed = 284651; % make sampling reproducible.
sampler.spec.outputChainSize = 10000; % Use a small chain size for illustration.
sampler.spec.parallelismNumThread = []; % Use these many parallel threads.
sampler.run ( @(x) pm.stats.dist.himmelblau.getLogUDF(x(1), x(2)) ...
            , 2 ...
            );

chain = sampler.readChain();
p = pm.vis.plot.Scatter(chain.df, "coly", "proposalAdaptation");
p.make("axes", {"yscale", "log"});
p.savefig("proposalAdaptation.png", "-m4");chain = sampler.readChain();

p = pm.vis.tile.Scatter(chain.df, "tileshape", [2, 1]);
p.make("coly", chain.sampleLogFuncColIndex + [1 : 2]);
p.savefig("trace.png", "-m4");
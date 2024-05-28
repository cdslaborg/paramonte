cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";
sampler.spec.proposalStart = [-5, 5];
sampler.spec.outputFileName = "himmelblau";
sampler.spec.randomSeed = 28457353; % make sampling reproducible.
sampler.spec.outputChainSize = 30000; % Use a small chain size for illustration.
sampler.spec.parallelismNumThread = []; % Use these many parallel threads.
sampler.spec.outputRestartFileFormat = "ascii";
sampler.run ( @(x) pm.stats.dist.himmelblau.getLogUDF(x(1), x(2)) ...
            , 2 ...
            );

chain = sampler.readChain();
chain = chain{1};
p = pm.vis.plot.Scatter(chain.df, "coly", "proposalAdaptation");
p.make("axes", {"yscale", "log"});
p.savefig("proposalAdaptation.png", "-m4");

p = pm.vis.plot.LineScatter(chain.df, "colx", chain.sampleLogFuncColIndex + 1, "coly", chain.sampleLogFuncColIndex + 2);
p.make("colc", "sampleLogFunc");
p.savefig("domain.png", "-m4");

p = pm.vis.tile.Line(chain.df, "tileshape", [2, 1]);
p.make("coly", chain.sampleLogFuncColIndex + [1 : 2], "colc", "sampleLogFunc");
p.savefig("trace.png", "-m4");

restart = sampler.readRestart();
restart = restart{1};

p = pm.vis.plot.Ellipse3(restart.proposalCov, restart.proposalMean, restart.uniqueStateVisitCount');
p.make("axes", {"zscale", "log"});
p.savefig("proposalAdaptation.png", "-m4");
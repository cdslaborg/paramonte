cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

sampler = pm.sampling.Paradram();
sampler.spec.outputFileName = "./out/run";
sampler.spec.outputChainSize = 3000;
sampler.run ( @(x) -sum((x - [-5; 5]).^2) ...
            , 2 ...
            );
chain = sampler.readChainMarkov();
chain = chain{1};
head(chain.df)

p = pm.vis.PlotScatter(chain.df, "coly", "proposalAdaptation");
p.make("axes", {"yscale", "log"});
p.savefig("readChainMarkov.proposalAdaptation.png", "-m3");

p = pm.vis.PlotLineScatter(chain.df, "colx", chain.sampleLogFuncColIndex + 1, "coly", chain.sampleLogFuncColIndex + 2);
p.make("colc", "sampleLogFunc");
p.savefig("readChainMarkov.domain.png", "-m3");

p = pm.vis.TileLine(chain.df, "tileshape", [2, 1]);
p.make("coly", chain.sampleLogFuncColIndex + [1 : 2], "colc", "sampleLogFunc");
p.savefig("readChainMarkov.traceplot.png", "-m3");
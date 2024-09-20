close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

ndim = 4;
sampler = pm.sampling.Paradram();
sampler.spec.outputRestartFileFormat = "ascii";
sampler.spec.outputChainSize = 10000;
sampler.spec.outputFileName = "./test";
sampler.spec.outputStatus = "retry";
sampler.run(@(x) -sum(x.^2), ndim);
restart = sampler.readRestart();
sample = sampler.readSample();
cols = sample{1}.slfc + [1 : ndim];
tx = pm.vis.Triplex ( pm.vis.SubplotScatter(sample{1}.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                    , pm.vis.SubplotHistogram(sample{1}.df, "colx", cols) ...
                    , pm.vis.SubplotContour(sample{1}.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                    );
tx.layout.reset();
tx.layout.cbarh.enabled = true;
tx.layout.cbarv.enabled = true;
tx.layout.cbarh.position = [];
tx.layout.cbarv.position = [];
%tx.layout.tiling.position = [.15, .15, nan, nan];
%tx.layout.tiling.tile.width = [];
tx.make();
pm.vis.figure.savefig("Triplex.full.png", "-m3");
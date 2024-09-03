cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";
sampler.spec.outputFileName = "mvn";
sampler.spec.randomSeed = 28457353; % make sampling reproducible.
sampler.spec.outputChainSize = 30000; % Use a small chain size for illustration.
sampler.spec.outputRestartFileFormat = "ascii";
sampler.silent = true;
sampler.run ( @(x) -sum((x - [-5; 5; 10; -10]) .^ 2) ...
            , 4 ...
            );

restart = sampler.readRestart();
restart = restart{1};

pv = pm.vis.PlotEllipse3(restart.proposalCov, restart.proposalMean, transpose(restart.uniqueStateVisitCount));
pv.make("axes", {"zscale", "log"}, "dimx", [1, 3], "dimy", [1, 3] + 1);
pv.savefig("PlotEllipse3.1.png", "-m3");

close all;
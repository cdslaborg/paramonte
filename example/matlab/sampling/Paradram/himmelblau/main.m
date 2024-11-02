cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Setup the sampler.
%%%%

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";
sampler.spec.proposalStart = [-5, 5];
sampler.spec.outputFileName = "himmelblau";
sampler.spec.randomSeed = 28457353; % make sampling reproducible.
sampler.spec.outputChainSize = 30000; % Use a small chain size for illustration.
sampler.spec.parallelismNumThread = []; % Set this to a positive number to request that many parallel threads for the sampling.
sampler.spec.outputRestartFileFormat = "ascii";
% Set ``mpiname`` to ``pm.lib.mpi.choice()`` or your choice of MPI
% library ("intel", "openmpi", "mpich", ...) for MPI-parallel applications.
sampler.mpiname = ''; %pm.lib.mpi.choice();
sampler.silent = ~isempty(sampler.mpiname);

%%%%
%%%% Run the sampler.
%%%%

if  true % set to false to only post-process an existing simulation in the current folder.
    sampler.run ( @(x) pm.stats.dist.himmelblau.getLogUDF(x(1), x(2)) ...
                , 2 ...
                );
    if  pm.array.len(sampler.mpiname) > 0
        % We do not want to do any post-processing in distributed MPI-parallel mode.
        % It would be at least redundant but more importantly, potentially troubling.
        return;
    end
end

%%%%
%%%% Postprocess the output chain file.
%%%%

chain = sampler.readChain();
chain = chain{1};

%%%%
%%%% Make plots of proposal adaptation.
%%%%

p = pm.vis.PlotScatter(chain.df, "coly", "proposalAdaptation");
p.make("axes", {"yscale", "log"});
p.savefig("Paradram.himmelblau.proposalAdaptation.png", "-m3");

chain.vis.proposalAdaptation.line.make()
chain.vis.proposalAdaptation.line.savefig("Paradram.himmelblau.proposalAdaptation.line.png", "-m3");

chain.vis.proposalAdaptation.scatter.make()
chain.vis.proposalAdaptation.scatter.savefig("Paradram.himmelblau.proposalAdaptation.scatter.png", "-m3");

%%%%
%%%% Make triplet (corner) plots from the chain file.
%%%%

chain.vis.triplex.lshc2.make();
p.savefig("Paradram.himmelblau.triplex.lshc2.png", "-m3");

chain.vis.triplex.lshc3.make();
p.savefig("Paradram.himmelblau.triplex.lshc3.png", "-m3");

chain.vis.triplex.lshcf.make();
p.savefig("Paradram.himmelblau.triplex.lshcf.png", "-m3");

%%%%
%%%% The number `chain.slfc` corresponds to the data column with header "sampleLogFunc"`.
%%%%

p = pm.vis.PlotLineScatter(chain.df, "colx", chain.slfc + 1, "coly", chain.slfc + 2);
p.make("colc", "sampleLogFunc");
p.savefig("Paradram.himmelblau.domain.2d.png", "-m3");

p = pm.vis.PlotScatter3(chain.df, "colx", chain.slfc + 1, "coly", chain.slfc + 2, "colz", chain.slfc, "colc", chain.slfc);
p.make();
p.savefig("Paradram.himmelblau.domain.3d.png", "-m3");

p = pm.vis.TileLine(chain.df, "tileshape", [2, 1]);
p.make("coly", chain.slfc + [1 : 2], "colc", "sampleLogFunc");
p.savefig("Paradram.himmelblau.traceplot.png", "-m3");

%%%%
%%%% Postprocess the output restart file.
%%%%

restart = sampler.readRestart();
restart = restart{1};

p = pm.vis.PlotEllipse3(restart.proposalCov, restart.proposalMean, transpose(restart.uniqueStateVisitCount));
p.make("axes", {"zscale", "log"});
p.savefig("Paradram.himmelblau.proposalCov.png", "-m3");

%%%%
%%%% Postprocess the output report file.
%%%%

report = sampler.readReport();
report = report{1};

for partype = ["optimal", "perfect"]
    report.stats.parallelism.(partype).scaling.strong.vis.lineScatter.make();
    report.stats.parallelism.(partype).scaling.strong.vis.lineScatter.savefig("Paradram.himmelblau.parallelism." + partype + ".scaling.strong.png", "-m3");
    %p = pm.vis.PlotLineScatter(report.stats.parallelism.(partype).scaling.strong.value, "colx", "1", "coly", "2", "colc", "2");
    %p.make("axes", {"xscale", "log"}, "plot", {"linewidth", 2}, "scatter", {"size", 7});
    %p.savefig("Paradram.himmelblau.parallelism." + partype + ".scaling.strong.png", "-m3");
end
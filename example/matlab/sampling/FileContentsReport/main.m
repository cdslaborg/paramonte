cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Setup the sampler.
%%%%

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";

%%%%
%%%% Run the sampler.
%%%%

sampler.run ( @(x) pm.stats.dist.himmelblau.getLogUDF(x(1), x(2)) ...
            , 2 ...
            );

%%%%
%%%% Ensure postprocessing the output report file is done by only
%%%% one parallel process if distributed (MPI) parallelism is enabled.
%%%%

if  pm.lib.mpi.runtime.rankp1() == 1

    report = sampler.readReport();
    report = report{1};

    for parcond = ["sameeff", "zeroeff"]
        report.stats.parallelism.speedup.scaling.strong.(parcond).vis.lineScatter.make();
        report.stats.parallelism.speedup.scaling.strong.(parcond).vis.lineScatter.savefig("Paradram.himmelblau.parallelism.speedup.scaling.strong." + parcond + ".png", "-m3");
    end

end
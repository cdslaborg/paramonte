cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Setup the sampler.
%%%%

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";
% Set ``mpiname`` to ``pm.lib.mpi.choice()`` or your choice of MPI
% library ("intel", "openmpi", "mpich", ...) for MPI-parallel applications.
sampler.mpiname = ''; %pm.lib.mpi.choice();

%%%%
%%%% Run the sampler.
%%%%

sampler.run ( @(x) pm.stats.dist.himmelblau.getLogUDF(x(1), x(2)) ...
            , 2 ...
            );
if  pm.array.len(sampler.mpiname) > 0
    % We do not want to do any post-processing in distributed MPI-parallel mode.
    % It would be at least redundant but more importantly, potentially troubling.
    return;
end

%%%%
%%%% Postprocess the output report file.
%%%%

report = sampler.readReport();
report = report{1};

for parcond = ["sameeff", "zeroeff"]
    report.stats.parallelism.speedup.scaling.strong.(parcond).vis.lineScatter.make();
    report.stats.parallelism.speedup.scaling.strong.(parcond).vis.lineScatter.savefig("Paradram.himmelblau.parallelism.speedup.scaling.strong." + parcond + ".png", "-m3");
end
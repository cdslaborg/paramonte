close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Setup the sampler.
%%%%

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";
sampler.spec.proposalStart = [-5, 5];
sampler.spec.outputFileName = "himmelblau";
sampler.spec.randomSeed = 28457353; % make sampling reproducible.
sampler.spec.outputRestartFileFormat = "ascii";

%%%%
%%%% Run the sampler.
%%%%

sampler.run ( @(x) pm.stats.dist.himmelblau.getLogUDF(x(1), x(2)) ...
            , 2 ... ndim
            );

%%%%
%%%% Postprocess the output restart file.
%%%%

restart = sampler.readRestart();
restart = restart{1};

disp("class(restart)");
disp( class(restart) );

for vistype = ["line"]
    restart.vis.cascade.(vistype).make();
    restart.vis.cascade.(vistype).savefig   ( "FileContentsRestartDRAM." ...
                                            + restart.vis.cascade.(vistype).template.subplot.coly ...
                                            + "." + vistype + ".png" ...
                                            , "-m3");
end

for field = ["proposalCor", "proposalCov"]
    for vistype = ["ellipse", "ellipse3"]
        restart.vis.(field).cascade.(vistype).make();
        restart.vis.(field).cascade.(vistype).savefig   ( "FileContentsRestartDRAM." + field + "." + vistype + "." ...
                                                        + string(restart.vis.(field).cascade.ellipse3.template.subplot.dimx) ...
                                                        + string(restart.vis.(field).cascade.ellipse3.template.subplot.dimy) ...
                                                        + ".png" ...
                                                        , "-m3" ...
                                                        );
    end
end
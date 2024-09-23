close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

url = "https://raw.githubusercontent.com/cdslaborg/paramontex/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/";
file = "regression_powerlaw_data_paradram/out/regression_powerlaw_process_1_sample.txt";
sampler = pm.sampling.Paradram();
sample = sampler.readSample(url + file);
sample = sample{1};

for vistype = string(fields(sample.vis))'
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    for plottype = string(fields(sample.vis.(vistype)))'
        sample.vis.(vistype).(plottype).make("figure", {"visible", "off"});
        figname = join(["FileContentsSample", vistype, plottype], ".");
        if  strcmpi(vistype, "cascade")
            figname = [figname + "." + string(1 : numel(sample.vis.(vistype).(plottype).window))];
        end
        disp("figname")
        disp( figname )
        sample.vis.(vistype).(plottype).savefig( figname + ".png", "-m3");
    end
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
end
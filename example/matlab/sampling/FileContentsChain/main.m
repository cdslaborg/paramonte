close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

url = "https://raw.githubusercontent.com/cdslaborg/paramontex/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/";
file = "regression_powerlaw_data_paradram/out/regression_powerlaw_process_1_chain.txt";
sampler = pm.sampling.Paradram();
chain = sampler.readSample(url + file);
chain = chain{1};

for vistype = string(fields(chain.vis))'
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    for plottype = string(fields(chain.vis.(vistype)))'
        chain.vis.(vistype).(plottype).make("figure", {"visible", "off"});
        figname = join(["FileContentsChain", vistype, plottype], ".");
        if  strcmpi(vistype, "cascade")
            figname = figname + "." + string(1 : numel(chain.vis.(vistype).(plottype).window));
        end
        disp("figname")
        disp( figname )
        chain.vis.(vistype).(plottype).savefig(figname + ".png", "-m3");
    end
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
end
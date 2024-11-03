close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

fname = "regression_powerlaw_process_1_chain.txt";
path = "regression_powerlaw_data_paradram/out/" + fname;
url = "https://raw.githubusercontent.com/cdslaborg/paramontex/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/";
chain = pm.sampling.FileContentsChainDRAM(websave(fname, url + path));

for visfield = string(fields(chain.vis))'

    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    plottypes = [];
    if  contains(lower(visfield), "adaptation")
        varname = "proposalAdaptation";
        plottypes = ["line", "scatter"];
    elseif contains(lower(visfield), "burnin")
        varname = "burninLocation";
        plottypes = ["line"];
    elseif contains(lower(visfield), "mean")
        varname = "meanAcceptanceRate";
        plottypes = ["line"];
    end

    if ~isempty(plottypes)
        for plottype = plottypes
            chain.vis.(visfield).(plottype).make();
            chain.vis.(visfield).(plottype).savefig(join(["FileContentsChainDRAM", varname, plottype], ".") + ".png", "-m3");
        end
    end

    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

end
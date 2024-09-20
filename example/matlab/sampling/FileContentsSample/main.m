close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

url = "https://raw.githubusercontent.com/cdslaborg/paramontex/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/";
file = "regression_powerlaw_data_paradram/out/regression_powerlaw_process_1_chain.txt";
sampler = pm.sampling.Paradram();
sample = sampler.readSample(url + file); 
sample = sample{1};

%for vistype = "tile"%string(fields(sample.vis))'
%    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
%    disp("vistype")
%    disp( vistype )
%    for plottype = string(fields(sample.vis.(vistype)))'
%        %clear sampler;
%        disp("plottype")
%        disp( plottype )
%        sample.vis.(vistype).(plottype).make(); %"figure", {"Visible", "off"}
%        uicontrol('Visible', 'off')
%        %sample.vis.(vistype).(plottype).savefig([], "-m3");
%        %sample.vis.(vistype).(plottype).reset();
%    end
%    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
%end
%%figfile = join(["FileContentsSample", vistype, plottype, "png"], ".");
%%disp("figfile")
%%disp( figfile )
%%sample.vis.(vistype).(plottype).savefigs(figfile, "-m3");

sample.vis.tile.line.make();
sample.vis.tile.scatter.make();
sample.vis.tile.lineScatter.make();

sample.vis.tile.line3.make();
sample.vis.tile.scatter3.make();
sample.vis.tile.lineScatter3.make();

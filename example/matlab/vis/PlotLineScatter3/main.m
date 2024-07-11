cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

pv = pm.vis.PlotLineScatter3(df);
pv.subplot.title.titletext = "LineScatter3 Plot";
pv.make("colx", 1, "coly", 2, "colz", 3);
pv.savefig("PlotLineScatter3.1.png", "-m3");

close all;
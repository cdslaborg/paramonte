cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

pv = pm.vis.PlotScatter(df);
pv.subplot.title.titletext = "Scatter Plot";
pv.make("colx", 1, "coly", 2, "colc", 3);
pv.savefig("PlotScatter.1.png", "-m3");

close all;
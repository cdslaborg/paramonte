cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

dfcor = pm.stats.Cor(df(:, 1 : 2));
pv = pm.vis.PlotHeatmap(dfcor.val);
pv.subplot.title.titletxt = "Heatmap Plot";
pv.make();
pv.subplot.setColorLim();
pv.savefig("PlotHeatmap.1.png", "-m3");

close all;
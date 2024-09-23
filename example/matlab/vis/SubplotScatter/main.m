cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

figure("color", "white");
sv = pm.vis.SubplotScatter(df);
sv.title.titletext = "Scatter Subplot";
sv.make("colx", 1, "coly", 2, "colc", 3);
pm.vis.figure.savefig("SubplotScatter.1.png", "-m3");

close all;
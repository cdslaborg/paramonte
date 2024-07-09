cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table([r .* cos(theta); r .* sin(theta); r(end) - r]');
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

figure("color", "white");
sv = pm.vis.SubplotLineScatter(df);
sv.title.titletext = "LineScatter Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("SubplotLineScatter.1.png", "-m3");

close all;
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

figure("color", "white");
sv = pm.vis.SubplotHistfit(df);
sv.title.titletext = "Histfit Subplot";
sv.make("colx", 1);
pm.vis.figure.savefig("SubplotHistfit.1.png", "-m3");

close all;
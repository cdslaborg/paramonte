cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table([r .* cos(theta); r .* sin(theta); r(end) - r]');
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

figure("color", "white");
dfcor = pm.stats.Cor(df(:, 1 : 2));
sv = pm.vis.subplot.Subplot("Heatmap", dfcor.val);
sv.make();
title("Heatmap Subplot");
pm.vis.figure.savefig("Subplot.Heatmap.1.png", "-m3");

close all;
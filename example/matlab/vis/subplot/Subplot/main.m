cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table([r .* cos(theta); r .* sin(theta); r(end) - r]');
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

figure("color", "white");
sv = pm.vis.subplot.Subplot("Line", df);
sv.title.titletext = "Line Subplot";
sv.surface.lineWidth = 2;
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("Subplot.Line.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Line3", df);
sv.title.titletext = "Line3 Subplot";
sv.surface.lineWidth = 2;
sv.make("colx", 1, "coly", 2, "colz", 3);
pm.vis.figure.savefig("Subplot.Line3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Scatter", df);
sv.title.titletext = "Scatter Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("Subplot.Scatter.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Scatter3", df);
sv.title.titletext = "Scatter3 Subplot";
sv.make("colx", 1, "coly", 2, "colz", 3);
pm.vis.figure.savefig("Subplot.Scatter3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("LineScatter", df);
sv.title.titletext = "LineScatter Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("Subplot.LineScatter.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("LineScatter3", df);
sv.title.titletext = "LineScatter3 Subplot";
sv.make("colx", 1, "coly", 2, "colz", 3);
pm.vis.figure.savefig("Subplot.LineScatter3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Histogram", df);
sv.title.titletext = "Histogram Subplot";
sv.make("colx", 1);
pm.vis.figure.savefig("Subplot.Histogram.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Histogram2", df);
sv.title.titletext = "Histogram2 Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("Subplot.Histogram2.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Histfit", df);
sv.title.titletext = "Histfit Subplot";
sv.make("colx", 1);
pm.vis.figure.savefig("Subplot.Histfit.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Contour", df);
sv.title.titletext = "Contour Subplot";
sv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pm.vis.figure.savefig("Subplot.Contour2.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Contour3", df);
sv.title.titletext = "Contour3 Subplot";
sv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pm.vis.figure.savefig("Subplot.Contour3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.subplot.Subplot("Contourf", df);
sv.title.titletext = "Contourf Subplot";
sv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pm.vis.figure.savefig("Subplot.Contour2.1.png", "-m3");

figure("color", "white");
dfcor = pm.stats.Cor(df(:, 1 : 2));
sv = pm.vis.subplot.Subplot("Heatmap", dfcor.val);
sv.make();
title("Heatmap Subplot");
pm.vis.figure.savefig("Subplot.Heatmap.1.png", "-m3");

close all;
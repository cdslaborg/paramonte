cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

figure("color", "white");
sv = pm.vis.Subplot("Line", df);
sv.title.titletext = "Line Subplot";
sv.surface.lineWidth = 2;
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("SubplotLine.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Line3", df);
sv.title.titletext = "Line3 Subplot";
sv.surface.lineWidth = 2;
sv.make("colx", 1, "coly", 2, "colz", 3);
pm.vis.figure.savefig("SubplotLine3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Scatter", df);
sv.title.titletext = "Scatter Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("SubplotScatter.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Scatter3", df);
sv.title.titletext = "Scatter3 Subplot";
sv.make("colx", 1, "coly", 2, "colz", 3);
pm.vis.figure.savefig("SubplotScatter3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("LineScatter", df);
sv.title.titletext = "LineScatter Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("SubplotLineScatter.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("LineScatter3", df);
sv.title.titletext = "LineScatter3 Subplot";
sv.make("colx", 1, "coly", 2, "colz", 3);
pm.vis.figure.savefig("SubplotLineScatter3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Histogram", df);
sv.title.titletext = "Histogram Subplot";
sv.make("colx", 1);
pm.vis.figure.savefig("SubplotHistogram.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Histogram2", df);
sv.title.titletext = "Histogram2 Subplot";
sv.make("colx", 1, "coly", 2);
pm.vis.figure.savefig("SubplotHistogram2.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Histfit", df);
sv.title.titletext = "Histfit Subplot";
sv.make("colx", 1);
pm.vis.figure.savefig("SubplotHistfit.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Contour", df);
sv.title.titletext = "Contour Subplot";
sv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pm.vis.figure.savefig("SubplotContour.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Contour3", df);
sv.title.titletext = "Contour3 Subplot";
sv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pm.vis.figure.savefig("SubplotContour3.1.png", "-m3");

figure("color", "white");
sv = pm.vis.Subplot("Contourf", df);
sv.title.titletext = "Contourf Subplot";
sv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pm.vis.figure.savefig("SubplotContourf.1.png", "-m3");

figure("color", "white");
dfcor = pm.stats.Cor(df(:, 1 : 2));
sv = pm.vis.Subplot("Heatmap", dfcor.val);
sv.make();
title("Heatmap Subplot");
pm.vis.figure.savefig("SubplotHeatmap.1.png", "-m3");

%%%% Make Ellipse plots.

sampler = pm.sampling.Paradram();
sampler.spec.outputStatus = "retry";
sampler.spec.outputFileName = "mvn";
sampler.spec.randomSeed = 28457353; % make sampling reproducible.
sampler.spec.outputChainSize = 30000; % Use a small chain size for illustration.
sampler.spec.outputRestartFileFormat = "ascii";
sampler.silent = true;
sampler.run(@(x) -sum((x - [-5; 5; 10; -10]) .^ 2), 4);
restart = sampler.readRestart();
restart = restart{1};

figure("color", "white");
sv = pm.vis.SubplotEllipse(restart.proposalCov, restart.proposalMean);
pv.subplot.title.titletext = "Ellipse Subplot";
sv.make("axes", {"zscale", "log"}, "dimx", [1, 3], "dimy", [1, 3] + 1);
pm.vis.figure.savefig("SubplotEllipse.1.png", "-m3");

figure("color", "white");
sv = pm.vis.SubplotEllipse3(restart.proposalCov, restart.proposalMean, transpose(restart.uniqueStateVisitCount));
pv.subplot.title.titletext = "Ellipse3 Subplot";
sv.make("axes", {"zscale", "log"}, "dimx", [1, 3], "dimy", [1, 3] + 1);
pm.vis.figure.savefig("SubplotEllipse3.1.png", "-m3");

close all;
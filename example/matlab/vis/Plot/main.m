cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

theta = linspace(0, 8 * pi, 500);
r = 1.5 * theta;
df = array2table(ctranspose([r .* cos(theta); r .* sin(theta); r(end) - r]));
df.Properties.VariableNames = ["X = r cos(theta)", "Y = r sin(theta)", "Z = r"];

pv = pm.vis.Plot(pm.vis.Subplot("Line", df));
pv.subplot.title.titletext = "Line Plot";
pv.subplot.surface.lineWidth = 2;
pv.make("colx", 1, "coly", 2);
pv.savefig("PlotLine.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Line3", df));
pv.subplot.title.titletext = "Line3 Plot";
pv.subplot.surface.lineWidth = 2;
pv.make("colx", 1, "coly", 2, "colz", 3);
pv.savefig("PlotLine3.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Scatter", df));
pv.subplot.title.titletext = "Scatter Plot";
pv.make("colx", 1, "coly", 2);
pv.savefig("PlotScatter.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Scatter3", df));
pv.subplot.title.titletext = "Scatter3 Plot";
pv.make("colx", 1, "coly", 2, "colz", 3);
pv.savefig("PlotScatter3.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("LineScatter", df));
pv.subplot.title.titletext = "LineScatter Plot";
pv.make("colx", 1, "coly", 2);
pv.savefig("PlotLineScatter.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("LineScatter3", df));
pv.subplot.title.titletext = "LineScatter3 Plot";
pv.make("colx", 1, "coly", 2, "colz", 3);
pv.savefig("PlotLineScatter3.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Histogram", df));
pv.subplot.title.titletext = "Histogram Plot";
pv.make("colx", 1);
pv.savefig("PlotHistogram.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Histogram2", df));
pv.subplot.title.titletext = "Histogram2 Plot";
pv.make("colx", 1, "coly", 2);
pv.savefig("PlotHistogram2.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Histfit", df));
pv.subplot.title.titletext = "Histfit Plot";
pv.make("colx", 1);
pv.savefig("PlotHistfit.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Contour", df));
pv.subplot.title.titletext = "Contour Plot";
pv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pv.savefig("PlotContour.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Contour3", df));
pv.subplot.title.titletext = "Contour3 Plot";
pv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pv.savefig("PlotContour3.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.Subplot("Contourf", df));
pv.subplot.title.titletext = "Contourf Plot";
pv.make("colx", 1, "coly", 2, "xlim", [-20, 20], "ylim", [-20, 20]);
pv.savefig("PlotContourf.1.png", "-m3");

dfcor = pm.stats.Cor(df(:, 1 : 2));
pv = pm.vis.Plot(pm.vis.Subplot("Heatmap", dfcor.val));
pv.make();
title("Heatmap Plot");
pv.savefig("PlotHeatmap.1.png", "-m3");

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

pv = pm.vis.Plot(pm.vis.SubplotEllipse(restart.proposalCov, restart.proposalMean));
pv.subplot.title.titletext = "Ellipse Plot";
pv.make("axes", {"zscale", "log"}, "dimx", [1, 3], "dimy", [1, 3] + 1);
pv.savefig("PlotEllipse.1.png", "-m3");

pv = pm.vis.Plot(pm.vis.SubplotEllipse3(restart.proposalCov, restart.proposalMean, transpose(restart.uniqueStateVisitCount)));
pv.subplot.title.titletext = "Ellipse3 Plot";
pv.make("axes", {"zscale", "log"}, "dimx", [1, 3], "dimy", [1, 3] + 1);
pv.savefig("PlotEllipse3.1.png", "-m3");

close all;
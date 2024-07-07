cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

df = pm.stats.dist.mvu.getRand(zeros(3, 1), chol(pm.stats.dist.cov.getRand(3), 'lower');, 5000)';
sv = pm.vis.subplot.Subplot(pm.vis.plot.Scatter(df, "colx", 1:2:5, "coly", 2:2:6));
sv.make();
sv.savefigs(["Subplot.window.1.png", "Subplot.window.2.png", "Subplot.window.3.png"], "-m3");

df = pm.stats.dist.mvu.getRand(zeros(3, 1), chol([1 .9 .9; .9 1 .9; .9 .9 1]), 10000)';
s = pm.vis.subplot.Subplot("scatter", df);
s.make("colx", 1, "coly", 2, "colc", 3);

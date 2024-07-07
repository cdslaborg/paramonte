cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

ndim = 10;
cv = pm.vis.cascade.Heatmap(pm.stats.dist.cov.getRand(ndim));
cv.make();
cv.savefigs(["Heatmap.window.1.png", "Heatmap.window.2.png", "Heatmap.window.3.png"], "-m3");
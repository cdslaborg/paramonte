cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

ndim = 10;
cv = pm.vis.CascadeHeatmap(pm.stats.dist.cov.getRand(ndim));
cv.make("rows", 1:ndim, "colx", {1:ndim, 1:2:ndim, 2:2:ndim}, "precision", 2);
cv.savefigs(["CascadeHeatmap.window.1.png", "CascadeHeatmap.window.2.png", "CascadeHeatmap.window.3.png"], "-m3");
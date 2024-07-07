cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

range = -6 : 0.01 : 6;
[x, y] = meshgrid(range, range);
z = exp(pm.stats.dist.himmelblau.getLogUDF(x, y, 1));
f = pm.vis.figure.Figure("figure", {"color", "white"});
f.make();
surf(x, y, z, log(z), "EdgeColor", "none");
set(gca, "zscale", "log");
zlabel("Z");
f.savefig("Figure.himmelblau.3d.png", "-m3");
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

ndim = 10;
mean = zeros(ndim, 1);
cholow = chol(pm.stats.dist.cov.getRand(ndim), 'lower');
df = ctranspose(pm.stats.dist.mvn.getRand(mean, cholow, 5000));
df = [df, pm.stats.dist.mvn.getLogPDF(ctranspose(df), mean, pm.matrix.inv(cholow))];
cv = pm.vis.CascadeLineScatter(df, "colx", 1:2:5, "coly", 2:2:6, "colc", ndim + 1);
cv.make();
cv.savefig(["CascadeLineScatter.window.1.png", "CascadeLineScatter.window.2.png", "CascadeLineScatter.window.3.png"], "-m3");
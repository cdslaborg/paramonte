cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

cholow = chol(pm.stats.dist.cov.getRand(10), 'lower');
df = transpose(pm.stats.dist.mvu.getRand(zeros(length(cholow), 1), cholow, 5000));
cv = pm.vis.CascadeHistfit(df, "colx", 1:2:5);
cv.make();
cv.savefig(["CascadeHistfit.window.1.png", "CascadeHistfit.window.2.png", "CascadeHistfit.window.3.png"], "-m3");
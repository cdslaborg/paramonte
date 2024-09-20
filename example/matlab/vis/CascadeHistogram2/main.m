cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

cholow = chol(pm.stats.dist.cov.getRand(10), 'lower');
df = transpose(pm.stats.dist.mvu.getRand(zeros(length(cholow), 1), cholow, 5000));
cv = pm.vis.CascadeHistogram2(df, "colx", 1:2:5, "coly", 2:2:6);
cv.make();
cv.savefig(["CascadeHistogram2.window.1.png", "CascadeHistogram2.window.2.png", "CascadeHistogram2.window.3.png"], "-m3");
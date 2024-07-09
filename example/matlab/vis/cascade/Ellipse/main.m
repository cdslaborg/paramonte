cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

ndim = 10;
nell = 100;
center = zeros(ndim, nell);
gramian = zeros(ndim, ndim, nell);
for iell = 1 : nell
    gramian(:, :, iell) = pm.stats.dist.cov.getRand(ndim, log(iell));
end
cv = pm.vis.CascadeEllipse(gramian, center, [], "dimx", 1:2:5, "dimy", 2:2:6);
cv.make();
cv.savefigs(["Ellipse.window.1.png", "Ellipse.window.2.png", "Ellipse.window.3.png"], "-m3");
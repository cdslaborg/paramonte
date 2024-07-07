cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

ndim = 10;
nell = 100;
center = zeros(ndim, nell);
gramian = zeros(ndim, ndim, nell);
for iell = 1 : nell
    gramian(:, :, iell) = pm.stats.dist.cov.getRand(ndim, log(iell));
end
cv = pm.vis.cascade.Ellipse3(gramian, center, [], [], "dims", [[1:2:5]', [2:2:6]']);
cv.make();
cv.savefigs(["Ellipse3.window.1.png", "Ellipse3.window.2.png", "Ellipse3.window.3.png"], "-m3");
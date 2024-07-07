cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('covmat = pm.stats.dist.cov.getRand(3)')
                covmat = pm.stats.dist.cov.getRand(3)
pm.matlab.show('cholow = chol(covmat, "lower")')
                cholow = chol(covmat, "lower")
pm.matlab.show('invmat = pm.matrix.inv(cholow)')
                invmat = pm.matrix.inv(cholow)
pm.matlab.show('invmat - inv(covmat)')
                invmat - inv(covmat)
assert(all(all(abs(invmat - inv(covmat)) < 1.e-12)));
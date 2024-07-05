cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('acf = pm.stats.AutoCorr(unifrnd(0, 1, 1000, 1));')
                acf = pm.stats.AutoCorr(unifrnd(0, 1, 1000, 1));
pm.matlab.show('p = pm.vis.plot.Line([acf.lag, acf.val.val]); p.make("colx", 1, "coly", 2); yline(acf.bnd); yline(-acf.bnd);')
                p = pm.vis.plot.Line([acf.lag, acf.val.val]); p.make("colx", 1, "coly", 2); yline(acf.bnd); yline(-acf.bnd);
pm.matlab.show('p.savefig("AutoCorr.unifrnd.png", "-m4");')
                p.savefig("AutoCorr.unifrnd.png", "-m4");
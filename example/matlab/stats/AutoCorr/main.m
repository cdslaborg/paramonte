cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('acf = pm.stats.AutoCorr(unifrnd(0, 1, 1000, 4));')
                acf = pm.stats.AutoCorr(unifrnd(0, 1, 1000, 4));

pm.matlab.show('acf.vis.plot.line.make(); %yline(acf.bnd); yline(-acf.bnd);')
                acf.vis.plot.line.make(); %yline(acf.bnd); yline(-acf.bnd);
pm.matlab.show('acf.vis.plot.line.savefig("AutoCorr.unifrnd.plot.line.png", "-m3");')
                acf.vis.plot.line.savefig("AutoCorr.unifrnd.plot.line.png", "-m3");

pm.matlab.show('acf.vis.tile.line.make();')
                acf.vis.tile.line.make();
pm.matlab.show('acf.vis.tile.line.savefig("AutoCorr.unifrnd.tile.line.png", "-m3");')
                acf.vis.tile.line.savefig("AutoCorr.unifrnd.tile.line.png", "-m3");

pm.matlab.show('acf.vis.cascade.line.make();')
                acf.vis.cascade.line.make();
pm.matlab.show('acf.vis.cascade.line.savefig("AutoCorr.unifrnd.cascade.line." + string([1, 2, 3, 4]) + ".png", "-m3");')
                acf.vis.cascade.line.savefig("AutoCorr.unifrnd.cascade.line." + string([1, 2, 3, 4]) + ".png", "-m3");

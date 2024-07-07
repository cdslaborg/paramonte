cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

% Make a positive-definite random matrix.
pm.matlab.show()
pm.matlab.show("cholow = chol(pm.stats.dist.cov.getRand(5), 'lower');")
                cholow = chol(pm.stats.dist.cov.getRand(5), 'lower');
pm.matlab.show("df = pm.stats.dist.mvn.getRand(zeros(length(cholow), 1), cholow, 5000)';")
                df = pm.stats.dist.mvn.getRand(zeros(length(cholow), 1), cholow, 5000)';

pm.matlab.show()
pm.matlab.show('c = pm.stats.Cov(df); c.val')
                c = pm.stats.Cov(df); c.val

pm.matlab.show()
pm.matlab.show('c = pm.stats.Cov(df, "pearson"); c.val')
                c = pm.stats.Cov(df, "pearson"); c.val

pm.matlab.show()
pm.matlab.show('c = pm.stats.Cov(df, "spearman"); c.val')
                c = pm.stats.Cov(df, "spearman"); c.val

pm.matlab.show()
pm.matlab.show('c = pm.stats.Cov(df, "kendall"); c.val')
                c = pm.stats.Cov(df, "kendall"); c.val

pm.matlab.show()
pm.matlab.show('p = pm.vis.plot.Heatmap(c.val); p.make("precision", 2); p.subplot.setColorLim();')
                p = pm.vis.plot.Heatmap(c.val); p.make("precision", 2); p.subplot.setColorLim();
pm.matlab.show('p.savefig("Cov.unifrnd.png", "-m3");')
                p.savefig("Cov.unifrnd.png", "-m3");
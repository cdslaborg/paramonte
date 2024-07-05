cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

% Make a positive-definite random matrix.
pm.matlab.show()
pm.matlab.show('randmat = diag(unifrnd(0, .9, 5, 1), 0);')
                randmat = diag(unifrnd(0, .9, 5, 1), 0);
pm.matlab.show('randmat(randmat==0) = unifrnd(-.9, .9, numel(randmat) - size(randmat, 1), 1);')
                randmat(randmat==0) = unifrnd(-.9, .9, numel(randmat) - size(randmat, 1), 1);
pm.matlab.show("covmat = tril(randmat)' * tril(randmat)")
                covmat = tril(randmat)' * tril(randmat)
pm.matlab.show("cholow = chol(covmat, 'lower');")
                cholow = chol(covmat, 'lower');
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
pm.matlab.show('p.savefig("Cov.unifrnd.png", "-m4");')
                p.savefig("Cov.unifrnd.png", "-m4");
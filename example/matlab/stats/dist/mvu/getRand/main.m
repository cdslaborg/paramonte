cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.stats.dist.mvu.getRand(zeros(2, 1))')
pm.matlab.show( pm.stats.dist.mvu.getRand(zeros(2, 1)) )

pm.matlab.show()
pm.matlab.show('figure("color", "white"); histogram(pm.stats.dist.mvu.getRand(zeros(2, 1), [], 10000));')
                figure("color", "white"); histogram(pm.stats.dist.mvu.getRand(zeros(2, 1), [], 10000));
pm.matlab.show('pm.vis.figure.savefig("mvu.getRand.hist1.png");')
                pm.vis.figure.savefig("mvu.getRand.hist1.png");

pm.matlab.show()
pm.matlab.show('figure("color", "white"); histogram(pm.stats.dist.mvu.getRand([-3, 3], chol([1 .5; .5, 1], "lower"), 10000));')
                figure("color", "white"); histogram(pm.stats.dist.mvu.getRand([-3, 3], chol([1 .5; .5, 1], "lower"), 10000));
pm.matlab.show('pm.vis.figure.savefig("mvu.getRand.hist2.png");')
                pm.vis.figure.savefig("mvu.getRand.hist2.png");

pm.matlab.show()
pm.matlab.show('rand = pm.stats.dist.mvu.getRand([-3, 3], chol([1 .9; .9, 1], "lower"), 10000);')
                rand = pm.stats.dist.mvu.getRand([-3, 3], chol([1 .9; .9, 1], "lower"), 10000);
pm.matlab.show('figure("color", "white"); scatter(rand(1, :), rand(2, :), 5, ".");')
                figure("color", "white"); scatter(rand(1, :), rand(2, :), 5, ".");
pm.matlab.show('pm.vis.figure.savefig("mvu.getRand.scatter1.png");')
                pm.vis.figure.savefig("mvu.getRand.scatter1.png");
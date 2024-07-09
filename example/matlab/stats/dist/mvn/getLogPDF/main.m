cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show("pdf = exp(pm.stats.dist.mvn.getLogPDF([1, 2, 3]'))")
                pdf = exp(pm.stats.dist.mvn.getLogPDF([1, 2, 3]'))
pm.matlab.show("pdf_ref = mvnpdf([1, 2, 3]')")
                pdf_ref = mvnpdf([1, 2, 3]')
assert(abs(pdf - pdf_ref) < 1.e-12)

pm.matlab.show()
pm.matlab.show('mean = [-3, 3]; cholow = chol([1 -.9; -.9, 1], "lower")')
                mean = [-3, 3]; cholow = chol([1 -.9; -.9, 1], "lower");
pm.matlab.show('sample = pm.stats.dist.mvn.getRand(mean, cholow, 1000);')
                sample = pm.stats.dist.mvn.getRand(mean, cholow, 1000);
pm.matlab.show("sample = array2table([sample', exp(pm.stats.dist.mvn.getLogPDF(sample, mean, pm.matrix.inv(cholow)))]);")
                sample = array2table([sample', exp(pm.stats.dist.mvn.getLogPDF(sample, mean, pm.matrix.inv(cholow)))]);
pm.matlab.show('sample.Properties.VariableNames = ["X", "Y", "MVN PDF"];')
                sample.Properties.VariableNames = ["X", "Y", "MVN PDF"];
pm.matlab.show('p = pm.vis.PlotScatter3(sample, "colx", 1, "coly", 2, "colz", 3, "colc", 3); p.make();')
                p = pm.vis.PlotScatter3(sample, "colx", 1, "coly", 2, "colz", 3, "colc", 3); p.make();
pm.matlab.show('pm.vis.figure.savefig("mvn.getLogPDF.scatter.3d.png", "-m4");')
                pm.vis.figure.savefig("mvn.getLogPDF.scatter.3d.png", "-m4");
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

x = exp(linspace(log(0.01), log(10), 1000))';

logPDF = zeros(numel(x), 6);
logPDF(:,1) = pm.stats.dist.gengamma.getLogPDF(x, 1.0, 0.5, 0.5);
logPDF(:,2) = pm.stats.dist.gengamma.getLogPDF(x, 2.0, 0.5, 1.0);
logPDF(:,3) = pm.stats.dist.gengamma.getLogPDF(x, 0.5, 2.0, 0.5);
logPDF(:,4) = pm.stats.dist.gengamma.getLogPDF(x, 0.2, 5.0, 0.2);
logPDF(:,5) = pm.stats.dist.gengamma.getLogPDF(x, .14, 7.0, .14);
logPDF(:,6) = pm.stats.dist.gengamma.getLogPDF(x, 2.0, 5.0, 0.3);

df = array2table([x, exp(logPDF)]);
p = pm.vis.PlotLine(df, "colx", 1, "coly", 2:length(df{1,:}), "plot", {"linewidth", 2});
p.subplot.colormap.enabled = false;
p.subplot.ylabel.txt = "PDF";
p.subplot.xlabel.txt = "X";
p.subplot.xlim = [0, 8.];
p.subplot.ylim = [0, .8];
p.make();
p.savefig("gengamma.getLogPDF.line.png", "-m4");
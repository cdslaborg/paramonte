cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

x = exp(linspace(log(0.01), log(10), 1000))';

cdf = zeros(numel(x), 6);
cdf(:,1) = pm.stats.dist.gengamma.getCDF(x, 1.0, 0.5, 0.5);
cdf(:,2) = pm.stats.dist.gengamma.getCDF(x, 2.0, 0.5, 1.0);
cdf(:,3) = pm.stats.dist.gengamma.getCDF(x, 0.5, 2.0, 0.5);
cdf(:,4) = pm.stats.dist.gengamma.getCDF(x, 0.2, 5.0, 0.2);
cdf(:,5) = pm.stats.dist.gengamma.getCDF(x, .14, 7.0, .14);
cdf(:,6) = pm.stats.dist.gengamma.getCDF(x, 2.0, 5.0, 0.3);

df = array2table([x, cdf]);
p = pm.vis.PlotLine(df, "colx", 1, "coly", 2:length(df{1,:}), "plot", {"linewidth", 2});
p.subplot.colormap.enabled = false;
p.subplot.ylabel.txt = "CDF";
p.subplot.xlabel.txt = "X";
p.make();
p.savefig("gengamma.getCDF.line.png", "-m4");
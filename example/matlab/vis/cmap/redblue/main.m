cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

figure("color", "white");
imagesc(peaks(500));
colormap(pm.vis.cmap.redblue(256));
colorbar;
pm.vis.figure.savefig("redblue.1.png", "-m3");

figure("color", "white");
load topo
imagesc(0:360, -90:90, topo); axis("xy");
colormap(pm.vis.cmap.redblue());
colorbar;
pm.vis.figure.savefig("redblue.2.png", "-m3");
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

figure("color", "white");
imagesc(peaks(500));
colormap(pm.vis.cmap.cold(256));
colorbar;
pm.vis.figure.savefig("cold.1.png", "-m3");

figure("color", "white");
load topo
imagesc(0:360, -90:90, topo); axis("xy");
colormap(pm.vis.cmap.cold());
colorbar;
pm.vis.figure.savefig("cold.2.png", "-m3");
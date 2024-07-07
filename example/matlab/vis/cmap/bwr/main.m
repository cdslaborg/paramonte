cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

figure("color", "white");
imagesc(peaks(500));
colormap(pm.vis.cmap.bwr(256));
colorbar;
pm.vis.figure.savefig("bwr.1.png", "-m3");

figure("color", "white");
imagesc(peaks(500), [0, 8]);
colormap(pm.vis.cmap.bwr());
colorbar;
pm.vis.figure.savefig("bwr.2.png", "-m3");

figure("color", "white");
imagesc(peaks(250), [-6, 0]);
colormap(pm.vis.cmap.bwr());
colorbar;
pm.vis.figure.savefig("bwr.3.png", "-m3");

figure("color", "white");
surf(peaks);
colormap(pm.vis.cmap.bwr());
axis("tight");
pm.vis.figure.savefig("bwr.4.png", "-m3");
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

figure("color", "white"); hold on; box on;
bcrd = pm.geom.ell2.getBorders();
h = plot(bcrd(:, 1), bcrd(:, 2), '-');
pm.vis.figure.savefig("getBorders.2d.png", "-m3");

npnt = 50;
figure("color", "white"); hold on; box on; view(3);
bcrd = pm.geom.ell2.getBorders([], [], repmat([1 : 20], npnt, 1));
for iell = 1 : 2 : size(bcrd, 2) / 3
    icol = (iell - 1) * 3 + 1;
    plot3(bcrd(:, icol), bcrd(:, icol + 1), bcrd(:, icol + 2), '-');
end
pm.vis.figure.savefig("getBorders.3d.png", "-m3");
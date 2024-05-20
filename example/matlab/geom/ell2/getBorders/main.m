cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

bcrd = pm.geom.ell2.getBorders();
figure; h = plot(bcrd(:, 1), bcrd(:, 2), '-');

bcrd = pm.geom.ell2.getBorders([], [], 20);
figure; h = plot3(bcrd(:, 1), bcrd(:, 2), bcrd(:, 3), '-');

npnt = 50;
figure; hold on; view(3);
bcrd = pm.geom.ell2.getBorders([], [], repmat([1 : 20], npnt, 1));
for iell = 1 : 2 : size(bcrd, 2) / 3
    icol = (iell - 1) * 3 + 1;
    plot3(bcrd(:, icol), bcrd(:, icol + 1), bcrd(:, icol + 2), '-');
end

%bcrd = pm.geom.ell2.getBorder();
%p = pm.vis.plot.Line(bcrd);
%p.make("colx", 1, "coly", 2);
%p.savefig("getBorder.2d.png", "-m4");
%
%bcrd = pm.geom.ell2.getBorder([], []);
%bcrd = [bcrd, zeros(size(bcrd, 1), 1)];
%p = pm.vis.plot.Line3(bcrd);
%p.make("colx", 1, "coly", 2, "colz", 3);
%p.savefig("getBorder.3d.png", "-m4");
%
%npnt = 500;
%range = 10 * pi * [-1 : 2 / (npnt - 1) : 1]';
%bcrd = [pm.geom.ell2.getBorder([], [], npnt), sin(range / 5) + sin(range) / 5];
%p = pm.vis.plot.Line3(bcrd);
%p.make("colx", 1, "coly", 2, "colz", 3);
%p.savefig("getBorder.wavy.png", "-m4");
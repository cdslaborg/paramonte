srcdir = fileparts(mfilename('fullpath'));
cd(srcdir); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

bcrd = pm.geom.ell2.getBorder();
p = pm.vis.plot.Line(bcrd);
p.make("colx", 1, "coly", 2, "figure", {"visible", "off", "color", "white"});
p.savefig("getBorder.2d.png", "-m4");

bcrd = pm.geom.ell2.getBorder([], []);
bcrd = [bcrd, zeros(size(bcrd, 1), 1)];
p = pm.vis.plot.Line3(bcrd);
p.make("colx", 1, "coly", 2, "colz", 3, "figure", {"visible", "off", "color", "white"});
p.savefig("getBorder.3d.png", "-m4");

npnt = 500;
range = 10 * pi * [-1 : 2 / (npnt - 1) : 1]';
bcrd = [pm.geom.ell2.getBorder([], [], npnt), sin(range / 5) + sin(range) / 5];
p = pm.vis.plot.Line3(bcrd);
p.make("colx", 1, "coly", 2, "colz", 3, "figure", {"visible", "off", "color", "white"});
p.savefig("getBorder.wavy.png", "-m4");
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show("figure('color', 'white');")
                figure('color', 'white');
pm.matlab.show("pm.vis.axes.getAxes(2, 1, 1, 'sv', 0, 'mr', 0); imagesc(magic(3));")
                pm.vis.axes.getAxes(2, 1, 1, 'sv', 0, 'mr', 0); imagesc(magic(3));
pm.matlab.show("pm.vis.axes.getAxes(2, 'pt', .02); imagesc(magic(4));")
                pm.vis.axes.getAxes(2, 'pt', .02); imagesc(magic(4));
pm.matlab.show("pm.vis.figure.savefig('getAxes.png', '-m4');")
                pm.vis.figure.savefig('getAxes.png', '-m4');
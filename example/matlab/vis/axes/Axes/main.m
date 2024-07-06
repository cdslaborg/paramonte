cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('a = pm.vis.axes.Axes("line")')
                a = pm.vis.axes.Axes("line")
pm.matlab.show()

pm.matlab.show()
pm.matlab.show('a = pm.vis.axes.Axes("line", "plot", {"linewidth", 2})')
                a = pm.vis.axes.Axes("line", "plot", {"linewidth", 2})
pm.matlab.show('hash = a.comp2hash("plot")')
                hash = a.comp2hash("plot")
pm.matlab.show()
pm.matlab.show('a.plot')
                a.plot
pm.matlab.show()
pm.matlab.show('a.reset("plot", {"markerSize", 5});')
                a.reset("plot", {"markerSize", 5});
pm.matlab.show()
pm.matlab.show('a.plot')
                a.plot
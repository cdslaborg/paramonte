cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show("val = repmat([5., 5.2, 5.5, 5.8, 6.], 20, 1);")
                val = repmat([5., 5.2, 5.5, 5.8, 6.], 20, 1);
pm.matlab.show("val")
pm.matlab.show( val )
pm.matlab.show("pm.math.pnint(val)")
pm.matlab.show( pm.math.pnint(val) )
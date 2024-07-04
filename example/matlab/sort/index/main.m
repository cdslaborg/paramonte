cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('array = [2.5, 4, 1.5, -3];')
                array = [2.5, 4, 1.5, -3];
pm.matlab.show('pm.sort.index(array)')
pm.matlab.show( pm.sort.index(array) )
pm.matlab.show('array(pm.sort.index(array))')
pm.matlab.show( array(pm.sort.index(array)) )

[~, indx] = sort(array);
assert(all(pm.sort.index(array) == indx));
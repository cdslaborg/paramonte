cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show("pm.math.logsubexp(log(400), log(1000))")
pm.matlab.show( pm.math.logsubexp(log(400), log(1000)) )

pm.matlab.show()
pm.matlab.show("pm.math.logsubexp(log(1000.), log(1001))")
pm.matlab.show( pm.math.logsubexp(log(1000.), log(1001)) )
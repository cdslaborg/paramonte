cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.web.isurl("https://github.com/cdslaborg/paramonte")')
                pm.web.isurl("https://github.com/cdslaborg/paramonte")

pm.matlab.show()
pm.matlab.show('pm.web.isurl("https://paramonte")')
                pm.web.isurl("https://paramonte")
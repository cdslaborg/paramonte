cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('hlink = pm.web.href("https://github.com/cdslaborg/paramonte")')
                hlink = pm.web.href("https://github.com/cdslaborg/paramonte")
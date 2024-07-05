cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('name = pm.web.basename("https://github.com/cdslaborg/paramonte")')
                name = pm.web.basename("https://github.com/cdslaborg/paramonte")

pm.matlab.show()
pm.matlab.show('name = pm.web.basename("https://github.com/cdslaborg/paramonte/")')
                name = pm.web.basename("https://github.com/cdslaborg/paramonte/")
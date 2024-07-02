cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('h = pm.matlab.Handle();')
                h = pm.matlab.Handle();

pm.matlab.show()
pm.matlab.show('h.help();')
                h.help();
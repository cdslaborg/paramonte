cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('weblinks = pm.lib.weblinks()')
                weblinks = pm.lib.weblinks();
pm.matlab.show("weblinks")
pm.matlab.show( weblinks )
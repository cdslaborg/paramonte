cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show('pm.matlab.isgui()')
pm.matlab.show( pm.matlab.isgui() )
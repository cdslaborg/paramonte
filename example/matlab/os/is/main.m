cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.os.is.lin()')
pm.matlab.show( pm.os.is.lin() )


pm.matlab.show()
pm.matlab.show('pm.os.is.mac()')
pm.matlab.show( pm.os.is.mac() )


pm.matlab.show()
pm.matlab.show('pm.os.is.win()')
pm.matlab.show( pm.os.is.win() )
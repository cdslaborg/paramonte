cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show("fid = fopen('temp.temp', 'w'); fclose(fid);")
                fid = fopen('temp.temp', 'w'); fclose(fid);

pm.matlab.show()
pm.matlab.show('pm.sys.path.list()')
pm.matlab.show( pm.sys.path.list() )

pm.matlab.show()
pm.matlab.show('pm.sys.path.rmfile("temp.temp")')
pm.matlab.show( pm.sys.path.rmfile("temp.temp") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.rmfile("temp.temp")')
pm.matlab.show( pm.sys.path.rmfile("temp.temp") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.rmfile("temp.temp", "Temp file is missing.")')
pm.matlab.show( pm.sys.path.rmfile("temp.temp", "Temp file is missing.") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.list()')
pm.matlab.show( pm.sys.path.list() )
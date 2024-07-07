cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.vis.color.rgb("amber")')
pm.matlab.show( pm.vis.color.rgb("amber") )

pm.matlab.show()
pm.matlab.show('pm.vis.color.rgb()')
pm.matlab.show( pm.vis.color.rgb() )
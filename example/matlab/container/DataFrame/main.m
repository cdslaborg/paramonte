cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

dfref = pm.container.DataFrame()
dfref.copy()

data = randi([0, 9], 10, 5)
dfref = pm.container.DataFrame(data)
dfref.copy()

dfref = pm.container.DataFrame(table(data))
dfref.copy()

dfref.ncol()
dfref.nrow()
dfref.rowslog()
dfref.rowslog(10)
dfref.rowslog(10, 5, 100)
dfref.rowslog(10, 5, 1000)
dfref.rowslog([], 5, 1000)
dfref.rowslog(10, -1, 1000)
dfref.rowslog(10, -1, -2)
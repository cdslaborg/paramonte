cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

df = pm.data.DataFrame()
df.copy()

data = randi([0, 9], 10, 5)
df = pm.data.DataFrame(data)
df.copy()

df = pm.data.DataFrame(table(data))
df.copy()

df.ncol()
df.nrow()
df.rowslog()
df.rowslog(10)
df.rowslog(10, 5, 100)
df.rowslog(10, 5, 1000)
df.rowslog([], 5, 1000)
df.rowslog(10, -1, 1000)
df.rowslog(10, -1, -2)
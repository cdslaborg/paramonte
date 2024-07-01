cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

df = pm.data.DataRef()
df.copy()

data = randi([0, 9], 10, 5)
df = pm.data.DataRef(data)
df.copy()

df = pm.data.DataRef(table(data))
df.copy()
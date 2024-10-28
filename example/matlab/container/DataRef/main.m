cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

df = pm.container.DataRef()
df.copy()

data = randi([0, 9], 10, 5)
df = pm.container.DataRef(data)
df.copy()

df = pm.container.DataRef(table(data))
df.copy()
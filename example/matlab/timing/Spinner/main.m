cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('spinner = pm.timing.Spinner();')
                spinner = pm.timing.Spinner();
pm.matlab.show('for i = 1 : 10; spinner.spin(i/10); end')
                for i = 1 : 10; spinner.spin(i/10); end
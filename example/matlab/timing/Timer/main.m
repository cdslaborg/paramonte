cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('timer = pm.timing.Timer();')
                timer = pm.timing.Timer();
pm.matlab.show('timer.tic(); % reset/start timer.')
                timer.tic(); % reset/start timer.
pm.matlab.show('pause(1);')
                pause(1);
pm.matlab.show('timer.toc()')
pm.matlab.show( timer.toc() )
pm.matlab.show('pause(.5);')
                pause(.5);
pm.matlab.show('timer.toc()')
pm.matlab.show( timer.toc() )
pm.matlab.show('pause(.2);')
                pause(.2);
pm.matlab.show('timer.del()')
pm.matlab.show( timer.del() )
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('linelist = split(pm.sys.bash_profile.getContents(), newline);')
                linelist = split(pm.sys.bash_profile.getContents(), newline);
pm.matlab.show('linelist(1:min(2, length(linelist)))') % Show only the first two lines for privacy reasons.
pm.matlab.show( linelist(1:min(2, length(linelist))) )
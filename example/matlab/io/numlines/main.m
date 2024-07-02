cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.io.numlines("main.m")')
disp( pm.io.numlines("main.m") )
disp("")


disp("")
disp('fid = fopen("main.m");')
      fid = fopen("main.m");
disp('line = fgetl(fid);')
      line = fgetl(fid);
disp('pm.io.numlines(fid)')
disp( pm.io.numlines(fid) )
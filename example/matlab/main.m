warning off;
srcpath = string(mfilename('fullpath')) + ".m";
rootdir = string(fileparts(srcpath));
cd(rootdir); % Change working directory to source code directory.
defpath = path; % startup path search list.

% Glob all example files to run.
if ~exist("examlist", "var")
    addpath('../'); % Add the ParaMonte library package to the MATLAB search list.
    examlist = pm.sys.path.glob(rootdir + "**main.m");
else
    % This is supposed to be passed to the routine 
    % through CMake call to ``matlab`` command.
    examlist = string(examlist);
    if  length(examlist) == 1
        examlist = strsplit(examlist, ";");
    end
end

% Loop over all examples to run.
for ipath = 1 : length(examlist)
    exampath = examlist(ipath);
    if ~strcmp(exampath, srcpath)
        cd(fileparts(exampath));
        outfile = "main.out.m";
        if isfile(outfile)
            delete(outfile);
        end
        diary main.out.m;
        disp("Running example: " + exampath);
        path(defpath);
        run(exampath);
        diary off;
        close all;
    end
end
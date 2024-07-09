warning off;
format compact;
srcpath = string(mfilename('fullpath')) + ".m";
rootdir = string(fileparts(srcpath));
cd(rootdir); % Change working directory to source code directory.
defpath = path; % startup path search list.

% Glob all example files to run.
if ~isfile("exam.list")
    addpath('../'); % Add the ParaMonte library package to the MATLAB search list.
    examlist = pm.sys.path.glob(rootdir + "**main.m");
else
    % The file is supposed to be created by CMake
    % or other tools (or human) and placed in the
    % root directory of the example folder.
    examdirs = readlines("exam.list");
    % Prepend the base directory to all directory paths.
    examlist = strings(length(examdirs), 1);
    for idir = 1 : length(examdirs)
        examlist(idir) = rootdir + "/" + string(examdirs(idir)) + "main.m";
        if ~isdir(examlist(idir))
            warning("The specified example directory does not exist: """ + examlist(idir) + """")
        end
    end
end

% Loop over all examples to run.
for ipath = 1 : length(examlist)
    exampath = examlist(ipath);
    if ~strcmp(exampath, srcpath)
        cd(fileparts(exampath));
        outfile = "main.out.m";
        if  isfile(outfile)
            delete(outfile);
        end
        disp("Running example: " + exampath);
        diary main.out.m;
        path(defpath);
        runExample(exampath);
        diary off;
        close all;
    end
end

cd(rootdir);

function runExample(exampath)
    % This function ensures all example variables remain in local scopes
    % and do not interfere with variables and scopes of other examples.
    run(exampath);
end
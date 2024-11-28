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
        examfile = fullfile(rootdir, string(examdirs(idir)), "main.m");
        if ~isfile(examfile)
            warning("The corresponding example file for the specified example directory does not exist: """ + examfile + """")
        else
            examlist(idir) = examfile;
        end
    end
end

% Loop over all examples to run.
outfile = "main.out.m";
for ipath = 1 : length(examlist)
    exampath = examlist(ipath);
    if ~strcmp(exampath, srcpath)
        cd(fileparts(exampath));
        if  isfile(outfile)
            delete(outfile);
        end
        disp("Running example: " + exampath);
        diary(outfile);
        path(defpath);
        runExample(exampath);
        diary off;
        close all;
        sanitize(outfile);
    end
end

cd(rootdir);

function runExample(exampath)
    % This function ensures all example variables remain in local scopes
    % and do not interfere with variables and scopes of other examples.
    run(exampath);
end

%%%% ensure Carriage Return, Backspace, and other special characters
%%%% are correctly interpreted or at least taken as new line characters.
function sanitize(outfile)
    try
        output = fileread(outfile);
        searchlist = [string(char(8)), string(char(13))]; % Backspace, CR
        for search = searchlist
            double = search + search;
            while contains(output, double)
                output = strrep(output, double, search);
            end
            output = strrep(output, search, newline);
        end
        fid = fopen(outfile, "w");
        fprintf(fid, "%s", output);
        fclose(fid);
    end
end
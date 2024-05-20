warning off;
srcpath = string(mfilename('fullpath')) + ".m";
rootdir = string(fileparts(srcpath));
cd(rootdir); % Change working directory to source code directory.

% Glob all example files to run.
if ~exist("examlist", "var")
    examlist = pm.sys.path.glob(rootdir + "**main.m");
end

% Loop over all examples to run.
for ipath = 1 : length(examlist)
    exampath = examlist(ipath);
    if ~strcmp(exampath, srcpath)
        disp("Running example: " + exampath);
        run(exampath);
        diary off;
        close all;
    end
end
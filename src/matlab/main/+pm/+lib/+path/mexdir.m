%
%   Return the vector of MATLAB strings containing the directory
%   path(s) containing the specified ParaMonte library mex shared file.
%
%       mexname
%
%           The input scalar MATLAB string containing
%           the (partial) name of a MEX file name without MEX extension.
%           (**optional**. If missing, all detected paths will be returned.)
%
%       config
%
%           The input scalar (or vector of) MATLAB string(s) containing the
%           (set of) keyword(s) to match in the identified MEX files paths.
%           (**optional**. If missing, all detected paths will be returned.)
%
%>  \return
%       mexdirs
%
%           The output vector of MATLAB strings containing the **unique**
%           directory path(s) to the ParaMonte library MEX shared files.
%           The output ``mexdirs`` will be an empty list
%           if no such paths are found.
%
%   Interface
%   ---------
%
%       mexdirs = pm.lib.path.mexdir()
%       mexdirs = pm.lib.path.mexdir([])
%       mexdirs = pm.lib.path.mexdir([], [])
%       mexdirs = pm.lib.path.mexdir(mexname)
%       mexdirs = pm.lib.path.mexdir([], config)
%       mexdirs = pm.lib.path.mexdir(mexname, [])
%       mexdirs = pm.lib.path.mexdir(mexname, config)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function mexdirs = mex(mexname, config)
    if  nargin < 2
        config = [];
    end
    if  nargin < 1
        mexname = [];
    end
    mexpaths = pm.lib.path.mex(mexname, config);
    mexdirs = strings(length(mexpaths), 1);
    for ipath = 1 : length(mexpaths)
        [dirpath, ~, ~] = fileparts(mexpaths{ipath});
        mexdirs(ipath) = string(dirpath);
    end
    mexdirs = unique(mexdirs);
end
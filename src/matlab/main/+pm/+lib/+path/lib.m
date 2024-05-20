%
%   Return a scalar MATLAB string containing the path to
%   the ``lib`` directory of the ParaMonte library package.
%
%       None
%
%>  \return
%       path
%
%           The output scalar MATLAB string containing the path to
%           the ``lib`` directory of the ParaMonte library package.
%
%   Interface
%   ---------
%
%       path = pm.lib.path.lib()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function path = lib(key)
    persistent path_persistent
    if isempty(path_persistent)
        path_persistent = string(fullfile(pm.lib.path.root(), "lib"));
    end
    path = path_persistent;
end
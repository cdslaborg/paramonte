%>  \brief
%>  Return a scalar MATLAB string containing the
%>  root directory of the ParaMonte library package.
%>
%>  \return
%>  ``path``    :   A scalar MATLAB string containing the root
%>                  directory of the ParaMonte library package.<br>
%>
%>  \interface{root}
%>  \code{.m}
%>
%>      path = pm.lib.path.root()
%>
%>  \endcode
%>
%>  \example{root}
%>  \include{lineno} example/lib/path/root/main.m
%>  \output{root}
%>  \include{lineno} example/lib/path/root/main.out.m
%>
%>  \final{root}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:51 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function path = root()
    persistent path_persistent
    if isempty(path_persistent)
        [path_persistent, ~, ~] = fileparts(mfilename('fullpath'));
        path_persistent = string(pm.sys.path.abs(fullfile(path_persistent, "..", "..")));
    end
    path = path_persistent;
    %addpath(genpath(tree.root),'-begin');
end
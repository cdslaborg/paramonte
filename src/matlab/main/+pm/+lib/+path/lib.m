%>  \brief
%>  Return a scalar MATLAB string containing the path to
%>  the ``lib`` directory of the ParaMonte library package.
%>
%>  \return
%>  ``path``    :   The output scalar MATLAB string containing the path to
%>                  the ``lib`` directory of the ParaMonte library package.
%>
%>  \interface{lib}
%>  \code{.m}
%>
%>      path = pm.lib.path.lib();
%>
%>  \endcode
%>
%>  \example{lib}
%>  \include{lineno} example/lib/path/lib/main.m
%>  \output{lib}
%>  \include{lineno} example/lib/path/lib/main.out.m
%>
%>  \final{lib}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:39 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function path = lib()
    persistent path_persistent
    if isempty(path_persistent)
        path_persistent = string(fullfile(pm.lib.path.root(), "lib"));
    end
    path = path_persistent;
end
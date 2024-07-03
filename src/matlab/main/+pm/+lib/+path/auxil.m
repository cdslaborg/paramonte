%>  \brief
%>  Return a scalar MATLAB string containing the path to
%>  the ``auxil`` directory of the ParaMonte library package.
%>
%>  \return
%>  ``path``    :   The output scalar MATLAB string containing the path to
%>                  the ``auxil`` directory of the ParaMonte library package.
%>
%>  \interface{auxil}
%>  \code{.m}
%>
%>      path = pm.lib.path.auxil()
%>
%>  \endcode
%>
%>  \example{auxil}
%>  \include{lineno} example/lib/path/auxil/main.m
%>  \output{auxil}
%>  \include{lineno} example/lib/path/auxil/main.out.m
%>
%>  \final{auxil}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:29 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function path = auxil()
    persistent path_persistent
    if isempty(path_persistent)
        path_persistent = string(fullfile(pm.lib.path.root(), "auxil"));
    end
    path = path_persistent;
end
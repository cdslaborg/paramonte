%>  \brief
%>  Return a MATLAB string containing the lower-case name of the current OS.
%>
%>  \return
%>  ``str`` :   The output MATLAB string that is either:<br>
%>              <ol>
%>                  <li>    ``"linux"`` if the OS is Linux.
%>                  <li>    ``"windows"`` if the OS is Windows.
%>                  <li>    ``"darwin"`` if the OS is macOS (Darwin).
%>              </ol>
%>
%>  \interface{namel}
%>  \code{.m}
%>
%>      str = pm.os.namel()
%>
%>  \endcode
%>
%>  \example{namel}
%>  \include{lineno} example/os/namel/main.m
%>  \output{namel}
%>  \include{lineno} example/os/namel/main.out.m
%>
%>  \final{namel}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:53 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function result = namel()
    if ispc
        result = "windows";
    elseif ismac
        result = "darwin";
    elseif isunix
        result = "linux";
    end
end
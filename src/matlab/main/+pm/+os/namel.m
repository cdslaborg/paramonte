%>  \brief
%>  Return a MATLAB string containing the lower-case name of the current OS.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `str`   :   The output MATLAB string containing either:
%>              ``"linux"`` if the OS is Linux.
%>              ``"windows"`` if the OS is Windows.
%>              ``"darwin"`` if the OS is macOS (Darwin).
%>
%>  \interface{namel}
%>  \code{.m}
%>
%>      str = pm.os.namel()
%>
%>  \endcode
%>  \final{namel}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:53 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
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
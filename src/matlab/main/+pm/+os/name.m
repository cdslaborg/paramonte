%>  \brief
%>  Return a MATLAB string containing the name of the current OS.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `str`   :   The output MATLAB string that is either:
%>              ``"Linux"`` if the OS is Linux.
%>              ``"Windows"`` if the OS is Windows.
%>              ``"Darwin"`` if the OS is macOS (Darwin).
%>
%>  \interface{name}
%>  \code{.m}
%>
%>      str = pm.os.name()
%>
%>  \endcode
%>  \final{name}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:52 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = name()
    if ispc
        str = "Windows";
    elseif ismac
        str = "Darwin";
    elseif isunix
        str = "Linux";
    end
end
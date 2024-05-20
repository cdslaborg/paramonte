%
%   Return a scalar MATLAB string containing ``4``
%   blank characters equivalent to a tab character.
%
%   This function primarily exists to bring
%   consistency to the ParaMonte library IO tasks.
%
%       None
%
%>  \return
%       str
%
%           The output scalar MATLAB string containing ``4``
%           blank characters equivalent to a tab character.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       pm.io.tab()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = tab()
    str = "    ";
end
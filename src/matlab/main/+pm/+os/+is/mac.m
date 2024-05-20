%
%   Return ``true`` if the current OS is MacOS (Darwin).
%
%       None
%
%>  \return
%       itis
%
%           The output MATLAB logical scalar value that is ``true`` 
%           if and only if the OS is MacOS (Darwin), otherwise ``false``.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       itis = pm.os.is.mac()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function result = mac()
    result = ismac;
end
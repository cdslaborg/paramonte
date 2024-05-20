%
%   Return ``true`` if the current OS is Windows.
%
%       None
%
%>  \return
%       itis
%
%           The output MATLAB logical scalar value that is ``true``
%           if and only if the OS is Windows, otherwise ``false``.
%
%   Interface
%   ---------
%
%       itis = pm.os.is.win()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = win()
    itis = ispc;
end
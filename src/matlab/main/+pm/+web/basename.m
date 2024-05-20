%
%   Return a scalar MATLAB string containing
%   the basename part of an input ``url``, defined as
%   the segment of the ``url`` after the last separator ``/``.
%
%       url
%
%           The input scalar MATLAB string containing a bare URL.
%
%>  \return
%       name
%
%           The output scalar MATLAB string containing
%           the basename part of an input ``url``, defined as
%           the segment of the ``url`` after the last separator ``/``.
%           If the input ``url`` ends with ``/``, then the output
%           basename is set to the default ``index.html``.
%
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       name = pm.web.basename(url)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function name = basename(url)
    name = strsplit(url, '/');
    name = name(end);
end

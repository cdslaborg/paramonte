%
%   Return a list of MATLAB strings containing the names of
%   all supported compiler suites (vendor names) used for
%   building the ParaMonte MATLAB shared libraries.
%
%   \devnote
%
%       The CS names within this function must be
%       regularly updated with the latest CS names
%       available in the ParaMonte installation guide.
%
%       None
%
%>  \return
%       typelist
%
%           A MATLAB string list containing:
%
%               ``["intel", "gnu"]``
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       typelist = pm.lib.clstypes()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function typelist = clstypes()
    typelist = ["intel", "gnu"];
end
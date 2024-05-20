%
%   Return a list of MATLAB strings containing the names of
%   all currently possible builds of the ParaMonte MATLAB
%   shared libraries.
%
%   \devnote
%
%       The build names within this function must be
%       regularly updated with the latest build names
%       available in the ParaMonte installation guide.
%
%       None
%
%>  \return
%       typelist
%
%           A MATLAB string list containing:
%
%               ``["native", "tuned", "ipo", "release", "testing", "debug"]``
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       typelist = pm.lib.bldtypes()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function typelist = bldtypes()
    typelist = ["native", "tuned", "ipo", "release", "testing", "debug"];
end
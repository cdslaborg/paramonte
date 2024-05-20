%
%   Return a list of MATLAB strings containing the lower-case names
%   of all OS platforms supported by the ParaMonte MATLAB library.
%
%       None
%
%>  \return
%       names
%
%           The output MATLAB string list containing:
%
%           ``["windows", "linux", "darwin"]``
%
%   Interface
%   ---------
%
%       names = pm.os.listl()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function names = listl()
    names = lower(pm.os.list());
end
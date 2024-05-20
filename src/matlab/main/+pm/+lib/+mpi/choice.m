%
%   Return the ParaMonte-preferred MPI library vendor name
%   as used in naming the ParaMonte MATLAB shared libraries.
%   Support for the output MPI library name by this routine is
%   guaranteed on the current platform.
%
%       None
%
%>  \return
%       name
%
%           The output scalar MATLAB string containing the
%           the ParaMonte-preferred MPI library vendor name as
%           used in naming the ParaMonte MATLAB shared libraries.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       name = pm.lib.mpi.choice()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function name = choice()
    if ispc() || pm.os.is.lin()
        name = "impi";
    elseif ismac()
        name = "openmpi";
    end
end
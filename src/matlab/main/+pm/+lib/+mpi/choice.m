%>  \brief
%>  Return the ParaMonte-preferred MPI library vendor name
%>  as used in naming the ParaMonte MATLAB shared libraries.
%>
%>  \details
%>  Support for the output MPI library name by this
%>  routine is guaranteed on the current platform.
%>
%>  \return
%>  `name`  :   The output scalar MATLAB string containing the
%>              the ParaMonte-preferred MPI library vendor name as
%>              used in naming the ParaMonte MATLAB shared libraries.
%>
%>  \interface{choice}
%>  \code{.m}
%>
%>      name = pm.lib.mpi.choice()
%>
%>  \endcode
%>
%>  \example{choice}
%>  \include{lineno} example/lib/mpi/choice/main.m
%>  \output{choice}
%>  \include{lineno} example/lib/mpi/choice/main.out.m
%>
%>  \final{choice}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:16 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function name = choice()
    if ispc() || pm.os.is.lin()
        name = "impi";
    elseif ismac()
        name = "openmpi";
    end
end
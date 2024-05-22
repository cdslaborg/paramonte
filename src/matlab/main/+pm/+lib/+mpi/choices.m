%>  \brief
%>  Return the ParaMonte-preferred MPI library vendor name(s)
%>  as used in naming the ParaMonte MATLAB shared libraries
%>  in the order of preference on the current platform.
%>
%>  \note
%>  Only the default (first) mpi library name is guaranteed to be
%>  supported in any pre-built distribution of the ParaMonte library.
%>
%>  \param  'None'
%>
%>  \return
%>  names   :   The output vector of MATLAB strings containing the
%>              the ParaMonte-preferred MPI library vendor names as
%>              used in naming the ParaMonte MATLAB shared libraries.
%>              in the default order of preference.
%>
%>  \interface{choices}
%>  \code{.m}
%>
%>      names = pm.lib.mpi.choices()
%>
%>  \endcode
%>  \final{choices}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:18 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function names = choices()
    if ispc()
        names = "impi";
    elseif pm.os.is.lin()
        names = ["impi", "mpich", "openmpi"];
    elseif ismac()
        names = ["openmpi", "mpich"];
    end
end
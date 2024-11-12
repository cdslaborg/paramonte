%>  \brief
%>  Return the MPI library vendor name from the input MPI name based
%>  on the internal ParaMonte shared libraries naming conventions.<br>
%>
%>  \param[in]  name    :   The output scalar MATLAB string containing the MPI
%>                          library name corresponding to the output MPI ``vendor``.<br>
%>                          If the input ``name`` is not supported, the default
%>                          string ``"mpi"`` is returned as the value of ``str``.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``"mpi"``       : The corresponding output ``vendor`` name would be unknown and set to ``"MPI"``.<br>
%>                              <li>    ``"impi"``      : The corresponding output ``vendor`` name would be set to ``"Intel"``.<br>
%>                              <li>    ``"mpich"``     : The corresponding output ``vendor`` name would be set to ``"MPICH"``.<br>
%>                              <li>    ``"openmpi"``   : The corresponding output ``vendor`` name would be set to ``"OpenMPI"``.<br>
%>                          </ol>
%>                          Note that **all values are case-insensitive**.<br>
%>
%>  \return
%>  ``str``             :   The output scalar MATLAB string, containing the MPI
%>                          library vendor name supported by the ParaMonte library.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``"MPI"``, representing an unknown MPI library vendor.
%>                              <li>    ``"Intel"``, representing the Intel MPI library vendor.
%>                              <li>    ``"MPICH"``, representing the MPICH MPI library vendor.
%>                              <li>    ``"OpenMPI"``, representing the OpenMPI library vendor.
%>                          </ol>
%>
%>  \interface{vendor}
%>  \code{.m}
%>
%>      str = pm.lib.mpi.vendor(name)
%>
%>  \endcode
%>
%>  \see
%>  [pm.lib.mpi.name()](@ref name)<br>
%>  [pm.lib.mpi.vendor()](@ref vendor)<br>
%>  [pm.sys.path.mpiexec.vendor()](@ref vendor)<br>
%>
%>  \example{vendor}
%>  \include{lineno} example/lib/mpi/vendor/main.m
%>  \output{vendor}
%>  \include{lineno} example/lib/mpi/vendor/main.out.m
%>
%>  \final{vendor}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:20 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = vendor(name)
    if  nargin < 1
        name = [];
    end
    if  isempty(name)
        help("pm.lib.mpi.vendor");
        disp("name");
        disp( name );
        error   ( newline ...
                + "The input argument ``name`` must be present and non-empty." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end
    nameLower = lower(string(name));
    for i = length(nameLower) : -1 : 1
        if  strcmp(nameLower(i), "openmpi") || strcmp(nameLower(i), "open-mpi") || strcmp(nameLower(i), "ompi")
            str(i) = "OpenMPI";
        elseif strcmp(nameLower(i), "mpich") || strcmp(nameLower(i), "mmpi")
            str(i) = "MPICH";
        elseif strcmp(nameLower(i), "intel") || strcmp(nameLower(i), "impi")
            str(i) = "Intel";
        else
            str(i) = "MPI";
        end
    end
end
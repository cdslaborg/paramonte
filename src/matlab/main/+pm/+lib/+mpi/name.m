%>  \brief
%>  Return the MPI library name as used in
%>  naming the ParaMonte MATLAB shared library directories.<br>
%>
%>  \param[in]  vendor  :   The input scalar MATLAB string, containing the MPI
%>                          library vendor supported by the ParaMonte library.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``"Intel"``, representing the Intel MPI library.
%>                              <li>    ``"MPICH"``, representing the MPICH MPI library.
%>                              <li>    ``"OpenMPI"``, representing the OpenMPI library.
%>                          </ol>
%>                          Note that **all values are case-insensitive**.<br>
%>
%>  \return
%>  ``str``             :   The output scalar MATLAB string containing the MPI
%>                          library name corresponding to the input MPI ``vendor``.<br>
%>                          If the input ``vendor`` is not supported, the default
%>                          string ``"mpi"`` is returned as the value of ``str``.<br>
%>
%>  \interface{name}
%>  \code{.m}
%>
%>      str = pm.lib.mpi.name(vendor)
%>
%>  \endcode
%>
%>  \see
%>  [pm.lib.mpi.name()](@ref name)<br>
%>  [pm.lib.mpi.vendor()](@ref vendor)<br>
%>  [pm.sys.path.mpiexec.vendor()](@ref vendor)<br>
%>
%>  \example{name}
%>  \include{lineno} example/lib/mpi/name/main.m
%>  \output{name}
%>  \include{lineno} example/lib/mpi/name/main.out.m
%>
%>  \final{name}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:20 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = name(vendor)
    if  nargin < 1
        vendor = [];
    end
    if  isempty(vendor)
        help("pm.lib.mpi.name");
        disp("vendor");
        disp( vendor );
        error   ( newline ...
                + "The input argument ``vendor`` must be present and non-empty." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end
    vendorLower = lower(string(vendor));
    for i = length(vendorLower) : -1 : 1
        if  strcmp(vendorLower(i), "openmpi") || strcmp(vendorLower(i), "open-mpi") || strcmp(vendorLower(i), "ompi")
            str(i) = "openmpi";
        elseif strcmp(vendorLower(i), "mpich") || strcmp(vendorLower(i), "mmpi")
            str(i) = "mpich";
        elseif strcmp(vendorLower(i), "intel") || strcmp(vendorLower(i), "impi")
            str(i) = "impi";
        else
            str(i) = "mpi";
        end
    end
end
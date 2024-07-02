%>  \brief
%>  Return the MPI library name as used in
%>  naming the ParaMonte MATLAB shared libraries.
%>
%>  \param[in]  vendor  :   The input scalar MATLAB string, containing the MPI
%>                          library vendor supported by the ParaMonte library.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``Intel``, representing the Intel MPI library.
%>                              <li>    ``MPICH``, representing the MPICH MPI library.
%>                              <li>    ``OpenMPI``, representing the OpenMPI library.
%>                          </ol>
%>                          Note that **all values are case-insensitive**.<br>
%>
%>  \return
%>  `str`               :   The output scalar MATLAB string containing the MPI
%>                          library name corresponding to the input MPI ``vendor``.<br>
%>                          If the input ``vendor`` is not supported, the default
%>                          string ``"mpi"`` is returned as the value ``str``.<br>
%>
%>  \interface{name}
%>  \code{.m}
%>
%>      str = pm.lib.mpi.name(vendor)
%>
%>  \endcode
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
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = name(vendor)
    vendorLower = lower(string(vendor));
    for i = length(vendorLower) : -1 : 1
        if strcmp(vendorLower(i), "openmpi") || strcmp(vendorLower(i), "open-mpi") || strcmp(vendorLower(i), "ompi")
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
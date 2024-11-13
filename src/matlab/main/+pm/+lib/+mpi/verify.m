%>  \brief
%>  Verify the existence of ParaMonte-compatible MPI
%>  library installations on the current system.<br>
%>
%>  \param[in]  vendor  :   The input scalar MATLAB string, containing the MPI
%>                          library vendor supported by the ParaMonte library.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``OpenMPI``, representing the OpenMPI library.<br>
%>                              <li>    ``MPICH``, representing the MPICH MPI library.<br>
%>                              <li>    ``Intel``, representing the Intel MPI library.<br>
%>                              <li>    ``any``, representing any available MPI library.<br>
%>                                      Specifying this value can force the function to
%>                                      return as soon as any ParaMonte-compatible MPI library
%>                                      installation is detected on the current system.<br>
%>                              <li>    ``all``, representing all available MPI library.<br>
%>                                      Specifying this value can lead to a comprehensive
%>                                      search for all available ParaMonte-compatible
%>                                      MPI installations on the system, which may
%>                                      time-consuming on some platforms.<br>
%>                                      This comprehensive system-level search
%>                                      only if the initial shallow search for
%>                                      the ParaMonte-compatible MPI library
%>                                      installations fails.<br>
%>                          </ol>
%>                          or any other informal name returned by the ParaMonte
%>                          MATLAB function [pm.lib.mpi.name()](@ref name).<br>
%>                          Note that **all values are case-insensitive**.<br>
%>                          (**optional, default = ``"any"``.)
%>
%>  \return
%>  ``failed``          :   The output scalar MATLAB ``logical`` that is ``true``
%>                          if and only if the MPI library verification fails.<br>
%>
%>  \interface{verify}
%>  \code{.m}
%>
%>      failed = pm.lib.mpi.verify()
%>      failed = pm.lib.mpi.verify(vendor)
%>
%>  \endcode
%>
%>  \example{verify}
%>  \include{lineno} example/lib/mpi/verify/main.m
%>  \output{verify}
%>  \include{lineno} example/lib/mpi/verify/main.out.m
%>
%>  \final{verify}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:18 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function failed = verify(vendor)

    failed = false;

    if  nargin < 1
        vendor = [];
    end

    if  isempty(vendor)
        vendor = "";
    end

    % isall = strcmpi(vendor, "all");
    % isany = strcmpi(vendor, "any");
    % if  isany || isall
    %     vendors = ["Intel", "MPICH", "OpenMPI"];
    % else
    %     vendors = string(vendor);
    % end
    %
    % mpiFound = false(numel(vendors), 1);
    % for iven = numel(vendors)
    %
    %     mpiVendor = vendors(iven);
    %     mpiVendorLower = lower(mpiVendor);
    %
    %     disp("Checking for the " + mpiVendor + " MPI library installations...")
    %
    %     [mpiname, ~, ~] = pm.lib.mpi.runtime.detect(mpiVendor);
    %
    %     mpiFound(iven) = ~strcmp(mpiname, "");
    %     if  mpiFound(iven)
    %         disp(pm.io.tab + "An " + mpiVendor + " MPI (runtime) library installation possibly exists on the system.")
    %     else
    %         disp(pm.io.tab + "None detected.")
    %     end
    %
    % end

    %%%%
    %%%% Perform a brute-force search for MPI installations.
    %%%%

    % if ~any(mpiFound)
        if  pm.sys.path.mpiexec.verify(vendor);
            failed = true;
        end
    % end

end
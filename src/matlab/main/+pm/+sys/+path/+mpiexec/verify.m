%>  \brief
%>  Verify the existence of ParaMonte-compatible
%>  MPI library installations on the current system
%>  using a brute-force system path search.<br>
%>
%>  \details
%>  Because of the intensive nature of the search,
%>  this algorithm may take some time to accomplish
%>  the task on some platforms with large file systems.<br>
%>
%>  The installation detection is performed by searching for
%>  available `mpiexec` MPI-launcher binaries for different
%>  MPI library installations.<br>
%>
%>  \param[in]  vendor  :   The input scalar or vector MATLAB string, containing the
%>                          MPI library vendor supported by the ParaMonte library.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``OpenMPI``, representing the OpenMPI library.<br>
%>                              <li>    ``MPICH``, representing the MPICH MPI library.<br>
%>                              <li>    ``Intel``, representing the Intel MPI library.<br>
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
%>                          Note that **all values are case-insensitive**.<br>
%>                          (**optional, default = ``"all"``.)
%>
%>  \return
%>  ``failed``          :   The output scalar MATLAB ``logical`` that is ``true``
%>                          if and only if the ``mpiexec`` verification fails.<br>
%>
%>  \interface{verify}
%>  \code{.m}
%>
%>      failed = pm.sys.path.mpiexec.verify()
%>      failed = pm.sys.path.mpiexec.verify(vendor)
%>
%>  \endcode
%>
%>  \warning
%>  This function has side-effects by printing messages on screen.<br>
%>
%>  \example{verify}
%>  \include{lineno} example/sys/path/mpiexec/verify/main.m
%>  \output{verify}
%>  \include{lineno} example/sys/path/mpiexec/verify/main.out.m
%>
%>  \final{verify}
%>
%>  \author
%>  \AmirShahmoradi, 11:51 PM Monday, November 11, 2024, Dallas, TX<br>
function failed = verify(vendor)

    failed = false;

    if  nargin < 1
        vendor = [];
    end

    if  isempty(vendor)
        vendor = "";
    end

    if  pm.array.len(vendor) == 0 || strcmpi(vendor, "any") || strcmpi(vendor, "all")
        mpiVendorList = ["Intel", "MPICH", "OpenMPI"];
    else
        mpiVendorList = [string(vendor)];
    end

    %%%%
    %%%% Perform a brute-force search for MPI installations.
    %%%%

    mpiexecFound = false;
    mpiVendorList = ["Intel", "MPICH", "OpenMPI"];
    for mpiVendor = mpiVendorList
        mpiVendorLower = lower(mpiVendor);
        disp("Checking for the " + mpiVendor + " MPI library installations. This may take a bit...")
        mpiPathList = pm.sys.path.mpiexec.find(mpiVendorLower);
        if  isempty(mpiPathList)
            disp(pm.io.tab() + "None detected.")
        else
            mpiexecFound = true;
            disp(pm.io.tab() + "The following " + mpiVendor + " ``mpiexec`` binaries were detected:")
            for ipath = 1 : length(mpiPathList)
                disp(pm.io.tab() + pm.io.tab() + mpiPathList(ipath))
            end
        end
    end

    if  mpiexecFound
        disp( newline ...
            + "At least one MPI runtime library installation from vendors " + newline ...
            + "compatible with the ParaMonte library exists on your system." + newline ...
            );
    else
        disp( newline ...
            + "No MPI runtime libraries were detected on your system." + newline ...
            + "Even if you believe there is any, it is invisible to runtime shell." + newline ...
            );
        failed = pm.lib.mpi.install();
    end

    if ~failed
        msgfile = fullfile(pm.lib.path.auxil(), ".mpi.run.msg");
        try
            fprintf('%s\n\n', strrep(fileread(msgfile), char(13), ''));
        catch me
            failed = true;
            weblinks = pm.lib.weblinks();
            warning ( newline ...
                    + string(me.identifier) + " : " + string(me.message) + newline ...
                    + "Failed to read the contents of the following file:" + newline ...
                    + newline ...
                    + pm.io.tab() + """" + msgfile + """" + newline ...
                    + newline ...
                    + "The ParaMonte MATLAB library integrity appears compromised." + newline ...
                    + "You can download the latest version of the library from: " + newline ...
                    + newline ...
                    + pm.io.tab() + pm.web.href(weblinks.github.releases.url) + newline ...
                    + newline ...
                    );
        end
    end

end
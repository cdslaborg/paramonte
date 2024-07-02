%>  \brief
%>  Check the requirements of the installed ParaMonte MATLAB library.
%>
%>  \details
%>  This function takes no input arguments and returns nothing.
%>
%>  \interface{verify}
%>  \code{.m}
%>
%>      pm.lib.verify();
%>
%>  \endcode
%>
%>  \example{verify}
%>  \include{lineno} example/lib/verify/main.m
%>  \output{verify}
%>  \include{lineno} example/lib/verify/main.out.m
%>
%>  \final{verify}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:11 PM, University of Texas at Arlington<br>
function verify()

    %if isunix
    %    setenv('PATH', ['/usr/local/bin:', getenv('PATH')]);
    %end

    % Require MATLAB >2017a.

    fprintf("Checking MATLAB version >2017b rquirement for compatibility with ParaMonte...")

    matlabVersion = version("-release");
    matlabVersion = str2double(matlabVersion(1:4));
    if matlabVersion < 2017
        warning ( newline ...
                + "MATLAB R" + string(matlabVersion) + "a or newer is required for proper ParaMonte functionality. " ...
                + "Please install a compatible version of MATLAB." ...
                + newline ...
                );
    else
        disp("PASSED")
    end

    weblinks = pm.lib.weblinks();

    % Ensure supported architecture.

    fprintf("Checking the system architecture compatibility with the ParaMonte library...")

    if ~strcmp(["amd64", "arm64"], pm.sys.arch())
        warning ( newline ...
                + "The ParaMonte library shared files are incompatible with the current system architecture:" + newline ...
                + newline ...
                + pm.io.tab + pm.sys.arch() + newline ...
                + newline ...
                + "All library functionalities that rely on the shared binary files will be unavailable." + newline ...
                + "However, all generic ParaMonte MATLAB functionalities are still available." + newline ...
                + "If you need the ParaMonte MATLAB offloading capabilities for this architecture," + newline ...
                + "let the library developers know by opening a GitHub issue at:" + newline ...
                + newline ...
                + pm.io.tab + pm.web.href(weblinks.github.issues.url) + newline ...
                + newline ...
                + "Alternatively, you can also build the library for your system" + newline ...
                + "by following the ParaMonte MATLAB build installations at:" + newline ...
                + newline ...
                + pm.io.tab + pm.web.href(weblinks.github.url) + newline ...
                + newline ...
                );
    else
        disp("PASSED")
    end

    % Check for MPI installations.

    mpiVendorList = ["Intel", "MPICH", "OpenMPI"];
    mpiexecFound = false;
    for vendor = mpiVendorList
        vendorLower = lower(vendor);
        disp("Checking the " + vendorLower + " MPI library installations. This may take a bit...")
        mpiPathList = pm.sys.path.mpiexec.find(vendorLower);
        if  isempty(mpiPathList)
            disp("None detected.")
        else
            mpiexecFound = true;
            disp(pm.io.tab + "The following " + vendor + " `mpiexec` binaries were detected:")
            for ipath = 1 : length(mpiPathList)
                disp(pm.io.tab + pm.io.tab + mpiPathList(ipath))
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
    end

    msgfile1 = fullfile(pm.lib.path.auxil(), ".mpi.install.msg");
    msgfile2 = fullfile(pm.lib.path.auxil(), ".mpi.run.msg");
    try
        fprintf('%s\n\n', fileread(msgfile1));
        fprintf('%s\n\n', fileread(msgfile2));
    catch me
        warning ( newline ...
                + string(me.identifier) + " : " + string(me.message) + newline ...
                + "Failed to read the contents of either or both of the following files:" + newline ...
                + newline ...
                + pm.io.tab + """" + msgfile1 + """" + newline ...
                + pm.io.tab + """" + msgfile2 + """" + newline ...
                + newline ...
                + "The ParaMonte MATLAB library integrty appears compromised." + newline ...
                + "You can download the latest version of the library from: " + newline ...
                + newline ...
                + pm.io.tab + pm.web.href(weblinks.github.release.url) + newline ...
                + newline ...
                );
    end

end
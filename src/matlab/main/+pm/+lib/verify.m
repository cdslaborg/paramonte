%>  \brief
%>  Verify the existence of the requirements of the installed ParaMonte MATLAB library.
%>
%>  \details
%>  This function takes no input arguments and returns nothing.
%>
%>  \return
%>  ``failed``  :   The output scalar MATLAB ``logical`` that is ``true`` if and
%>                  only if the ParaMonte library components verification fails.<br>
%>
%>  \interface{verify}
%>  \code{.m}
%>
%>      failed = pm.lib.verify();
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
function failed = verify()

    failed = false;

    %if isunix
    %    setenv('PATH', ['/usr/local/bin:', getenv('PATH')]);
    %end

    %%%%
    %%%% Require MATLAB >2017a.
    %%%%

    matlabMinimum = "2023a";
    fprintf("Checking MATLAB version >" + matlabMinimum + " requirement for compatibility with ParaMonte...")

    matlabVersion = version("-release");
    matlabVersion = str2double(matlabVersion(1:4));
    if  2023 <= matlabVersion
        disp("PASSED")
    else
        failed = true;
        warning ( newline ...
                + "MATLAB R" + string(matlabMinimum) + "a or newer may be required for " ...
                + "proper functionality of some of the ParaMonte MATLAB library routines, " ...
                + "particularly, the ParaMonte samplers and other routines which rely compiled MEX files. " ...
                + "Please install a compatible version of MATLAB. " ...
                + "Proceeding with the current installation..." ...
                + newline ...
                );
    end

    weblinks = pm.lib.weblinks();

    %%%%
    %%%% Ensure supported architecture.
    %%%%

    fprintf("Checking the system architecture compatibility with the ParaMonte library...")

    if  any(strcmp(["amd64", "arm64"], pm.sys.arch()))
        disp("PASSED")
    else
        failed = true;
        warning ( newline ...
                + "The ParaMonte library shared files are incompatible with the current system architecture:" + newline ...
                + newline ...
                + pm.io.tab() + pm.sys.arch() + newline ...
                + newline ...
                + "All library functionalities that rely on the shared binary files will be unavailable." + newline ...
                + "However, all generic ParaMonte MATLAB functionalities are still available." + newline ...
                + "If you need the ParaMonte MATLAB offloading capabilities for this architecture," + newline ...
                + "let the library developers know by opening a GitHub issue at:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(weblinks.github.issues.url) + newline ...
                + newline ...
                + "Alternatively, you can also build the library for your system" + newline ...
                + "by following the ParaMonte MATLAB build installations at:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(weblinks.github.url) + newline ...
                + newline ...
                );
    end

    %%%%
    %%%% Check for MPI installations.
    %%%%

    if  pm.lib.mpi.verify()
        failed = true;
    end

end
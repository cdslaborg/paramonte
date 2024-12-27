%>  \brief
%>  Install or provide instructions to install the MPI
%>  library from the requested input MPI ``vendor``.<br>
%>
%>  \details
%>  Currently, this function can only install the
%>  Intel MPI runtime libraries on Linux or Windows.<br>
%>  For other MPI vendors, it can only provide installation instructions.<br>
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
%>                          (**optional, default = ``"Intel"`` on platforms supported by the Intel MPI library, otherwise ``"all"``.)
%>  \param[in]  isyes   :   The input scalar MATLAB ``logical`.<br>
%>                          If provided, it will serve as the answer to the question displayed on the prompt by
%>                          the function to confirm whether the user wants to install the MPI libraries or not.<br>
%>                          This functionality becomes important when interaction with users through the prompt is impossible.<br>
%>                          (**optional. If missing or empty, the user will be prompted for an answer.)
%>
%>  \return
%>  ``failed``          :   The output scalar MATLAB ``logical`` that is ``true`` if and only if
%>                          the MPI library installation or its instruction display fails.<br>
%>
%>  \interface{install}
%>  \code{.m}
%>
%>      failed = pm.lib.mpi.install();
%>      failed = pm.lib.mpi.install(vendor);
%>      failed = pm.lib.mpi.install(vendor, []);
%>      failed = pm.lib.mpi.install(vendor, isyes);
%>      failed = pm.lib.mpi.install([], isyes);
%>
%>  \endcode
%>
%>  \warning
%>  This function has side-effects by printing messages on screen.<br>
%>
%>  \example{install}
%>  \include{lineno} example/lib/mpi/install/main.m
%>  \output{install}
%>  \include{lineno} example/lib/mpi/install/main.out.m
%>
%>  \final{install}
%>
%>  \author
%>  \AmirShahmoradi, 11:51 PM Monday, November 11, 2024, Dallas, TX<br>
function failed = install(vendor, isyes)

    failed = false;

    if  nargin < 2
        isyes = [];
    end

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

    weblinks = pm.lib.weblinks();

    %%%%
    %%%% Read the MPI installation instructions.
    %%%%

    msgfile = fullfile(pm.lib.path.auxil(), ".mpi.install.msg");
    try
        fprintf('%s\n\n', strrep(fileread(msgfile), char(13), ''));
    catch me
        failed = true;
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

    %%%%
    %%%% Install the Intel runtime library if requested.
    %%%%

    if  any(contains(pm.lib.mpi.vendor(mpiVendorList), "Intel")) && ~ismac

        promptmsg = "Would you like to install the Intel MPI runtime libraries on this system now (requires internet access)? (y/n)";
        if ~isempty(isyes)
            disp(promptmsg);
            itis = isyes;
        else
            itis = pm.io.isyes(promptmsg);
        end

        if  itis

            auxpath = pm.lib.path.auxil()
            auxfile = fullfile(auxpath, ".deps." + pm.os.namel() + ".mpi.intel");

            try

                mpifile = strrep(strrep(fileread(auxfile), char(13), ''), char(10), '');
                binlink = weblinks.github.releases.download.auxil.url + "/" + mpifile;
                binpath = websave(fullfile(auxpath, mpifile), binlink);

                disp("The binary file is stored at: """ + binpath + """");

                disp("Executing binary file...");
                disp("Please use the default recommended installation settings if and when prompted.");

                if  ispc
                    [failed, cmdout] = system(binpath, "-echo");
                else
                    [failed, cmdout] = system("chmod +x " + binpath + " && " + binpath, "-echo");
                end
                failed = failed ~= 0;

                if  failed
                    warning ( newline ...
                            + "Failed to execute the Intel MPI installation binary file." + newline ...
                            + "Here is the terminal output:" + newline ...
                            + newline ...
                            + strjoin(string(cmdout), newline) + newline ...
                            + newline ...
                            + "You can install the MPI library by manually executing the binary located at:" + newline ...
                            + newline ...
                            + pm.io.tab() + binpath + newline ...
                            + newline ...
                            + "Skipping the Intel MPI installation..." + newline ...
                            + newline ...
                            );
                end

            catch me

                failed = true;
                warning ( newline ...
                        + string(me.identifier) + " : " + string(me.message) + newline ...
                        + "Failed to read the contents of the following file:" + newline ...
                        + newline ...
                        + pm.io.tab() + """" + auxfile + """" + newline ...
                        + newline ...
                        + "The ParaMonte MATLAB library integrity appears compromised." + newline ...
                        + "You can download the latest version of the library from: " + newline ...
                        + newline ...
                        + pm.io.tab() + pm.web.href(weblinks.github.releases.tag.auxil.url) + newline ...
                        + newline ...
                        );

            end

        end

    end

end
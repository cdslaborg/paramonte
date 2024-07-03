%>  \brief
%>  Return a scalar MATLAB string containing the download
%>  weblink to latest ParaMonte MATLAB library version
%>  or return an empty string ``""`` is there is none.
%>
%>  \param[in]  silent  :   The input scalar MATLAB logical.<br>
%>                          If ``true``, all descriptive messages on
%>                          the MATLAB command line will be suppressed.<br>
%>                          (**optional**, default = ``false``)
%>
%>  \return
%>  ``str``             :   The output scalar MATLAB string containing the weblink to the latest
%>                          available version of the ParaMonte library that is newer than the
%>                          existing ParaMonte library version on the current system.<br>
%>                          The output ``str`` will be set to empty string ``""``
%>                          if there is no newer version or the function fails.
%>
%>  \interface{weblink}
%>  \code{.m}
%>
%>      str = pm.lib.update.weblink()
%>
%>  \endcode
%>
%>  \example{weblink}
%>  \include{lineno} example/lib/update/weblink/main.m
%>  \output{weblink}
%>  \include{lineno} example/lib/update/weblink/main.out.m
%>
%>  \final{weblink}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:54 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = weblink(silent)
    if  nargin < 1
        silent = [];
    end
    if  isempty(silent)
        silent = false;
    end
    str = "";
    weblinks = pm.lib.weblinks();
    currentVersionString = pm.lib.version();
    latestVersionString = pm.lib.update.version(silent);
    if latestVersionString == ""
        if ~silent
            warning ( newline ...
                    + "Failed to capture the latest ParaMonte MATLAB library version the project repository." + newline ...
                    + "If you believe this behavior is unexpected, please report this issue at: " + newline ...
                    + newline ...
                    + pm.io.tab + pm.web.href(weblinks.github.issues.url) + newline ...
                    + newline ...
                    );
        end
        return
    else
        currentVersionTriplet = pm.str.split2real(currentVersionString);
        latestVersionTriplet = pm.str.split2real(latestVersionString);
        if all(currentVersionTriplet == latestVersionTriplet)
            if ~silent
                disp( newline ...
                    + "It appears that you have the latest version of the ParaMonte MATLAB library." + newline ...
                    + newline ...
                    );
            end
        else
            updateAvailable =   (latestVersionTriplet(1)  > currentVersionTriplet(1)) ...
                            ||  (latestVersionTriplet(1) == currentVersionTriplet(1) && latestVersionTriplet(2)  > currentVersionTriplet(2)) ...
                            ||  (latestVersionTriplet(1) == currentVersionTriplet(1) && latestVersionTriplet(2) == currentVersionTriplet(2) && latestVersionTriplet(3) > currentVersionTriplet(3));
            if updateAvailable
                str = weblinks.github.release.latest.url + "/libparamonte_matlab_" + pm.os.namel();
                if ispc()
                    str = str + ".zip";
                else
                    str = str + ".tar.gz";
                end
                if ~silent
                    disp( newline ...
                        + "A newer ParaMonte MATLAB library version (" + latestVersionString + ") is " + newline ...
                        + "available on GitHub or is in preparation for release. " + newline ...
                        + "The currently-installed version on your system is: " + currentVersionString + newline ...
                        + "You can download the latest version of the ParaMonte MATLAB library from " + newline ...
                        + newline ...
                        + pm.io.tab + pm.web.href(str) + newline ...
                        + newline ...
                        );
                end
            else
                if ~silent
                    disp( newline ...
                        + "Looks like you have a version of the ParaMonte MATLAB library (" + currentVersionString + ") " + newline ...
                        + "that is newer than the most recent release version (" + latestVersionString + ")." + newline ...
                        + "Get in touch with us at:" + newline ...
                        + newline ...
                        + pm.io.tab + pm.web.href(weblinks.github.url) + newline ...
                        + newline ...
                        + "We need your superpowers!" + newline ...
                        );
                end
            end
        end
    end
end
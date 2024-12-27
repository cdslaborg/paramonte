%>  \brief
%>  Return a scalar MATLAB string containing the latest
%>  available ParaMonte MATLAB version newer than the
%>  existing version on the current system.
%>
%>  \param[in]  silent  :   The input scalar MATLAB logical.<br>
%>                          If ``true``, all descriptive messages on
%>                          the MATLAB command line will be suppressed.<br>
%>                          (**optional**, default = ``false``)
%>
%>  \return
%>  ``str``             :   The output scalar MATLAB string containing the
%>                          latest available ParaMonte MATLAB version newer
%>                          than the existing version on the current system.<br>
%>                          The output ``str`` will be set to empty string ``""``
%>                          if there is no newer version or the function fails.<br>
%>
%>  \interface{version}
%>  \code{.m}
%>
%>      str = pm.lib.update.version()
%>      str = pm.lib.update.version(silent)
%>
%>  \endcode
%>
%>  \example{version}
%>  \include{lineno} example/lib/update/version/main.m
%>  \output{version}
%>  \include{lineno} example/lib/update/version/main.out.m
%>
%>  \final{version}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:57 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = version(silent)
    if  nargin < 1
        silent = [];
    end
    if  isempty(silent)
        silent = false;
    end
    versionLink = "https://raw.githubusercontent.com/cdslaborg/paramonte/main/src/matlab/VERSION.md";
    try
        versionFileLineList = strsplit(webread(versionLink), newline);
        str = string(versionFileLineList{1});
    catch me
        str = "";
        if ~silent
            weblinks = pm.lib.weblinks();
            warning ( newline ...
                    + string(me.identifier) + " : " + string(me.message) + newline ...
                    + "Failed to fetch the latest version number from the weblink:" + newline ...
                    + newline ...
                    + pm.io.tab() + pm.web.href(versionLink) + newline ...
                    + newline ...
                    + "Ensure MATLAB has access to the internet." + newline ...
                    + "Otherwise, the structure of the ParaMonte library project on GitHub might have changed." + newline ...
                    + "The current ParaMonte library on your system is " + pm.lib.version() + newline ...
                    + "You can check for any newer package releases on the ParaMonte GitHub repository release page:" + newline ...
                    + newline ...
                    + pm.io.tab() + pm.web.href(weblinks.github.releases.url) + newline ...
                    + newline ...
                    );
        end
    end
end
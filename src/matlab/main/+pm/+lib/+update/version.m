%
%   Return a scalar MATLAB string containing the latest
%   available ParaMonte MATLAB version newer than the
%   existing version on the current system.
%
%       silent
%
%           The input scalar MATLAB logical.
%           If ``true``, all descriptive messages on
%           the MATLAB command line will be suppressed.
%           (**optional**, default = ``false``)
%
%>  \return
%       str
%
%           The output scalar MATLAB string containing the
%           latest available ParaMonte MATLAB version newer
%           than the existing version on the current system.
%           The output ``str`` will be set to empty string ``""``
%           if there is no newer version or the function fails.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       str = pm.lib.update.version()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = version(silent)
    if  nargin < 1
        silent = [];
    end
    if  isempty(silent)
        silent = false;
    end
    weblinks = pm.lib.weblinks();
    versionLink = "https://raw.githubusercontent.com/cdslaborg/paramonte/main/src/matlab/VERSION.md";
    try
        versionFileLineList = strsplit(webread(versionLink), newline);
        str = string(versionFileLineList{1});
    catch
        str = "";
        if ~silent
            warning ( newline ...
                    + "Failed to fetch the latest version number from the weblink:" + newline ...
                    + newline ...
                    + pm.io.tab + pm.web.href(versionLink) + newline ...
                    + newline ...
                    + "Ensure MATLAB has access to the internet." + newline ...
                    + "Otherwise, the structure of the ParaMonte library project on GitHub might have changed." + newline ...
                    + "The current ParaMonte library on your system is " + pm.lib.version() + newline ...
                    + "You can check for any newer package releases on the ParaMonte GitHub repository release page:" + newline ...
                    + newline ...
                    + pm.io.tab + pm.web.href(weblinks.github.release.url) + newline ...
                    + newline ...
                    );
        end
    end
end
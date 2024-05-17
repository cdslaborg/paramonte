function str = version(silent)
    %
    %   Return a scalar MATLAB string containing the latest
    %   available ParaMonte MATLAB version newer than the
    %   existing version on the current system.
    %
    %   Parameters
    %   ----------
    %
    %       silent
    %
    %           The input scalar MATLAB logical.
    %           If ``true``, all descriptive messages on
    %           the MATLAB command line will be suppressed.
    %           (**optional**, default = ``false``)
    %
    %   Returns
    %   -------
    %
    %       str
    %
    %           The output scalar MATLAB string containing the
    %           latest available ParaMonte MATLAB version newer
    %           than the existing version on the current system.
    %           The output ``str`` will be set to empty string ``""``
    %           if there is no newer version or the function fails.
    %
    %   Interface
    %   ---------
    %
    %       str = pm.lib.update.version()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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
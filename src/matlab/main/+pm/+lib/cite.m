function str = cite(format)
    %
    %   Return a scalar MATLAB raw string or HTML-style hyperlink
    %   containing the weblink to ParaMonte library publications.
    %
    %   Parameters
    %   ----------
    %
    %       format
    %
    %           The input scalar MATLAB string indicating
    %           the format of the output we address.
    %
    %               1.  An input value of ``"raw"`` will return a raw web address.
    %               2.  An input value of ``"html"`` will return an HTML style web address.
    %
    %           (**optional**, default = ``"html"``)
    %
    %   Returns
    %   -------
    %
    %       str
    %
    %           The output scalar MATLAB string containing either,
    %
    %               -   a raw weblink to the ParaMonte library publications
    %                   if the library is used outside the MATLAB GUI interface, or,
    %
    %               -   an HTML-style hyperlink to the ParaMonte library publications
    %                   if the library is used within the MATLAB GUI interface.
    %
    %   Interface
    %   ---------
    %
    %       str = pm.lib.cite()
    %       str = pm.lib.cite(format)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    citefile = fullfile(pm.lib.path.auxil(), ".cite.link");
    try
        if  nargin < 1
            format = "html";
        end
        str = string(fileread(citefile));
        if strcmpi(format, "html")
            str = pm.web.href(str);
        end
    catch
        str = "";
        weblinks = pm.lib.weblinks();
        warning ( newline ...
                + "Failed to read the contents of the citation file:" + newline ...
                + newline ...
                + pm.io.tab + """" + citefile + """" + newline ...
                + newline ...
                + "The structure of the ParaMonte library appears compromised." + newline ...
                + "You can always a fresh latest version of the library from:" + newline ...
                + newline ...
                + pm.io.tab + pm.web.href(weblinks.github.release.url()) + newline ...
                + newline ...
                );
    end
end
function hlink = href(url)
    %
    %   Return an HTML-style decoration of the input URL
    %   if the ParaMonte MATLAB library is used in GUI,
    %   otherwise, return the input URL as is.
    %
    %   This functionality is important for properly displaying
    %   hyperlinks within the MATLAB command-prompt while avoiding
    %   the unnecessary HTML syntax clutter where it is not recognized.
    %
    %   Parameters
    %   ----------
    %
    %       url
    %
    %           The input scalar MATLAB string containing a bare url.
    %
    %   Returns
    %   -------
    %
    %       hlink
    %
    %           The output scalar MATLAB string containing the
    %           hyperlink (HTML-style-decorated) version of the input ``url``.
    %
    %   Interface
    %   ---------
    %
    %       hlink = pm.web.href()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    persistent guiEnabled
    if isempty(guiEnabled)
        guiEnabled = pm.matlab.isgui();
    end
    if guiEnabled
        hlink = "<a href=""" + url + """>" + url + "</a>";
    else
        hlink = url;
    end
end

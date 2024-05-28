%>  \brief
%>  Return a scalar MATLAB raw string or HTML-style hyperlink
%>  containing the weblink to ParaMonte library publications.
%>
%>  \param[in]  format  :   The input scalar MATLAB string indicating
%>                          the format of the output we address.
%>                              1.  An input value of ``"raw"`` will return a raw web address.
%>                              2.  An input value of ``"html"`` will return an HTML style web address.
%>                          (**optional**, default = ``"html"``)
%>
%>  \return
%>  `str`               :   The output scalar MATLAB string containing either,
%>                          -   a raw weblink to the ParaMonte library publications
%>                              if the library is used outside the MATLAB GUI interface, or,
%>
%>                          -   an HTML-style hyperlink to the ParaMonte library publications
%>                              if the library is used within the MATLAB GUI interface.
%>
%>  \interface{cite}
%>  \code{.m}
%>
%>      str = pm.lib.cite()
%>      str = pm.lib.cite(format)
%>
%>  \endcode
%>  \final{cite}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:05 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = cite(format)
    citefile = fullfile(pm.lib.path.auxil(), ".cite.link");
    try
        if  nargin < 1
            format = "html";
        end
        str = string(fileread(citefile));
        if strcmpi(format, "html")
            str = pm.web.href(str);
        end
    catch me
        str = "";
        weblinks = pm.lib.weblinks();
        warning ( newline ...
                + string(me.identifier) + " : " + string(me.message) + newline ...
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
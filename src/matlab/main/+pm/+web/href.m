%>  \brief
%>  Return an HTML-style decoration of the input URL
%>  if the ParaMonte MATLAB library is used in GUI,
%>  otherwise, return the input URL as is.
%>
%>  \details
%>  This functionality is important for properly displaying
%>  hyperlinks within the MATLAB command-prompt while avoiding
%>  the unnecessary HTML syntax clutter where it is not recognized.
%>
%>  \param[in]  url :   The input scalar MATLAB string containing a bare url.
%>
%>  \return
%>  `hlink`         :   The output scalar MATLAB string containing the
%>                      hyperlink (HTML-style-decorated) version of the input ``url``.
%>
%>  \interface{href}
%>  \code{.m}
%>
%>      hlink = pm.web.href()
%>
%>  \endcode
%>
%>  \final{href}
%>
%>  \author
%>  \JoshuaOsborne, May 22 2024, 7:49 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hlink = href(url)
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

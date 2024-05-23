%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the input URL string exists in the world wide web.
%>
%>  \param[in]  url :   The input scalar MATLAB string
%>                      whose existence as a URL is to be tested.
%>
%>  \return
%>  `itis`          :   The output scalar MATLAB logical that is ``true`` if and
%>                      only if the input URL string exists in the world wide web.
%>
%>  \interface{isurl}
%>  \code{.m}
%>
%>      itis = pm.web.isurl(url)
%>
%>  \endcode
%>
%>  \final{isurl}
%>
%>  \author
%>  \JoshuaOsborne, May 22 2024, 7:50 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = isurl(url)
    itis = false;
    try
        urlobj = java.net.URL(url); %create the url object
        % Get the proxy information using the MATLAB proxy API.
        proxy = com.mathworks.webproxy.WebproxyFactory.findProxyForURL(urlobj);
        % Open a connection to the urlobj.
        if isempty(proxy)
            urlConnection = urlobj.openConnection;
        else
            urlConnection = urlobj.openConnection(proxy);
        end
        % Try to start the input stream.
        inputStream = urlConnection.getInputStream;
        itis = true;
    end
end
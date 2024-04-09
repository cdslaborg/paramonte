function itis = isurl(url)
    %
    %   Return a scalar MATLAB logical that is ``true`` if and
    %   only if the input URL string exists in the world wide web.
    %
    %   Parameters
    %   ----------
    %
    %       url
    %
    %           The input scalar MATLAB string
    %           whose existence as a URL is to be tested.
    %
    %   Returns
    %   -------
    %
    %       itis
    %
    %           The output scalar MATLAB logical that is ``true`` if and
    %           only if the input URL string exists in the world wide web.
    %
    %   Interface
    %   ---------
    %
    %       itis = pm.web.isurl(url)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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
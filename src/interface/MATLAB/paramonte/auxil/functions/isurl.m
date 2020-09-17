function tf = urlexist(url)
    URL = java.net.URL(url); %create the URL object
    % Get the proxy information using the MATLAB proxy API.
    proxy = com.mathworks.webproxy.WebproxyFactory.findProxyForURL(URL); 
    % Open a connection to the URL.
    if isempty(proxy)
        urlConnection = URL.openConnection;
    else
        urlConnection = URL.openConnection(proxy);
    end
    % Try to start the input stream
    try
        inputStream = urlConnection.getInputStream;
        tf = true;
    catch
        tf = false;
    end
end
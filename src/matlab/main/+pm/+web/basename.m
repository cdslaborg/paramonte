function name = basename(url)
    %
    %   Return a scalar MATLAB string containing
    %   the basename part of an input ``url``, defined as
    %   the segment of the ``url`` after the last separator ``/``.
    %
    %   Parameters
    %   ----------
    %
    %       url
    %
    %           The input scalar MATLAB string containing a bare URL.
    %
    %   Returns
    %   -------
    %
    %       name
    %
    %           The output scalar MATLAB string containing
    %           the basename part of an input ``url``, defined as
    %           the segment of the ``url`` after the last separator ``/``.
    %           If the input ``url`` ends with ``/``, then the output
    %           basename is set to the default ``index.html``.
    %
    %
    %   Interface
    %   ---------
    %
    %       name = pm.web.basename(url)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    name = strsplit(url, '/');
    name = name(end);
end

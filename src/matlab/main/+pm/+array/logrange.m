function array = logrange(start, stop, sizemax)
    %
    %   Return a set of maximum ``sizemax`` unique integer spacings almost
    %   linearly spaced in the natural logarithmic scale between the 
    %   specified ``start`` and ``stop`` of the output range.
    %
    %   Parameters
    %   ----------
    %
    %       start
    %
    %           The input scalar MATLAB whole-number (integer)
    %           representing the starting point of the output range.
    %
    %       stop
    %
    %           The input scalar MATLAB whole-number (integer)
    %           representing the stopping point of the output range.
    %
    %       sizemax
    %
    %           The input scalar MATLAB whole-number (integer)
    %           representing the maximum size of the output range.
    %           Due to rounding operation involved in creating the
    %           output range, it is impossible to prespecify the
    %           output range size, only the maximum.
    %           (**optional**, default = ``100``)
    %
    %   Returns
    %   -------
    %
    %       array
    %
    %           The output vector of MATLAB real values containing
    %           the set of naturally logarithmically-spaced integer 
    %           values in the specified input range.
    %
    %   Interface
    %   ---------
    %
    %       array = pm.array.logrange(start, stop)
    %       array = pm.array.logrange(start, stop, [])
    %       array = pm.array.logrange(start, stop, sizemax)
    %
    %   Example
    %   -------
    %
    %       array = pm.array.logrange(1, 1000) % 10 11 12 13 14 15 16 17 18 19 20
    %       array = pm.array.logrange(1, 1000, 10) % 1 2 4 8 16 32 63 126 251 501 1000
    %       array = pm.array.logrange(1, 1000, 20) % 10  15
    %       array = pm.array.logrange(1, 1000, [])
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if  nargin < 3
        sizemax = [];
    end
    if  isempty(sizemax)
        sizemax = 100;
    end
    logub = log(stop);
    loglb = log(start);
    logskip = (logub - loglb) / sizemax;
    array = unique(round(pm.array.logspace(loglb, logub, logskip)));
end
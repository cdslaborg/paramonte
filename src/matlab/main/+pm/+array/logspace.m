%
%   Return a set of unique real spacings linearly-spaced
%   in the logarithmic scale in the input given base,
%   between the specified lower and upper bounds.
%
%       loglb
%
%           The input scalar MATLAB real containing
%           the natural logarithm of the lower bound
%           of the output logarithmically-linear spaced vector.
%
%       logub
%
%           The input scalar MATLAB real containing
%           the natural logarithm of the upper bound
%           of the output logarithmically-linear spaced vector.
%
%       logskip
%
%           The input scalar MATLAB real of value larger than ``1``
%           containing the natural logarithm of the spacing between
%           the natural logarithm of the output values.
%           (**optional**, default = ``(logub - loglb) / 100``)
%
%       base
%
%           The input scalar MATLAB real
%           containing the base of the logarithmic space.
%           (**optional**, default = ``exp(1)``)
%
%>  \return
%       array
%
%           The output vector of MATLAB real values containing
%           the set of logarithmically-spaced values in the
%           specified input range with the specified ``base``.
%
%   Interface
%   ---------
%
%       array = pm.array.logspace(loglb, logub)
%       array = pm.array.logspace(loglb, logub, logskip)
%       array = pm.array.logspace(loglb, logub, [], base)
%       array = pm.array.logspace(loglb, logub, logskip, base)
%
%   Example
%   -------
%
%       array = pm.array.logspace(log(10), log(20))
%       array = pm.array.logspace(log(10), log(20), log(1.5)) % 10.000000000000002  15.000000000000007
%       array = pm.array.logspace(log(10), log(20), [], 2)
%       array = pm.array.logspace(log(10), log(20), log(1.5), 3) % 10.000000000000004  15.000000000000012
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function array = logspace(loglb, logub, logskip, base)
    if nargin < 4
        base = [];
    end
    if nargin < 3
        logskip = [];
    end
    if isempty(logskip)
        logskip = (logub - loglb) / 100;
    end
    if ~isempty(base)
        array = base .^ ((loglb : logskip : logub));
        %array = base .^ ((loglb : logskip : logub) ./ log(base));
    else
        array = exp(loglb : logskip : logub);
    end
end
%
%   Return ``log(exp(larger) - exp(smaller))`` more accurately.
%
%   \warning
%
%       The onus is on the user to ensure ``smaller < larger``.
%
%       smaller
%
%           The input scalar MATLAB real
%           representing the natural logarithm
%           of the smaller value.
%
%       smaller
%
%           The input scalar MATLAB real
%           representing the natural logarithm
%           of the larger value.
%
%>  \return
%       val
%
%           The output scalar MATLAB real containing the
%           the result of ``log(exp(larger) - exp(smaller))`` accurately.
%
%   Interface
%   ---------
%
%       val = pm.math.logsubexp(smaller, larger)
%
%   Example
%   -------
%
%       pm.math.logsubexp(log(400), log(1000))
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function val = logsubexp(smaller, larger)
    val = larger + log(-expm1(smaller - larger));
end
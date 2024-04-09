function val = logsubexp(smaller, larger)
    %
    %   Return ``log(exp(larger) - exp(smaller))`` more accurately.
    %
    %   \warning
    %
    %       The onus is on the user to ensure ``smaller < larger``.
    %
    %   Parameters
    %   ----------
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
    %   Returns
    %   -------
    %
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
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    val = larger + log(-expm1(smaller - larger));
end
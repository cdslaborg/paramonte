%>  \brief
%>  Return a set of unique real spacings linearly-spaced
%>  in the logarithmic scale in the input given base,
%>  between the specified lower and upper bounds.
%>
%>  \param[in]  loglb   :   The input scalar MATLAB real containing
%>                          the natural logarithm of the lower bound
%>                          of the output logarithmically-linear spaced vector.
%>  \param[in]  logub   :   The input scalar MATLAB real containing
%>                          the natural logarithm of the upper bound
%>                          of the output logarithmically-linear spaced vector.
%>  \param[in]  logskip :   The input scalar MATLAB real of value larger than ``1``
%>                          containing the natural logarithm of the spacing between
%>                          the natural logarithm of the output values.<br>
%>                          (**optional**, default = ``(logub - loglb) / 100``)
%>  \param[in]  base    :   The input scalar MATLAB real
%>                          containing the base of the logarithmic space.<br>
%>                          (**optional**, default = ``exp(1)``)
%>
%>  \return
%>  ``array``           :   The output vector of MATLAB real values containing
%>                          the set of logarithmically-spaced values in the
%>                          specified input range with the specified ``base``.
%>
%>  \interface{logspace}
%>  \code{.m}
%>
%>      array = pm.array.logspace(loglb, logub)
%>      array = pm.array.logspace(loglb, logub, logskip)
%>      array = pm.array.logspace(loglb, logub, [], base)
%>      array = pm.array.logspace(loglb, logub, logskip, base)
%>
%>  \endcode
%>
%>  \example{logspace}
%>  \include{lineno} example/array/logspace/main.m
%>  \output{logspace}
%>  \include{lineno} example/array/logspace/main.out.m
%>
%>  \final{logspace}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:33 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
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
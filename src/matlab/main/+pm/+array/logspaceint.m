%>  \brief
%>  Return a set of unique integer spacings almost linearly
%>  spaced in the logarithmic scale in the input given base,
%>  between the specified lower and upper bounds.
%>
%>  \param[in]  loglb   :   See the documentation of the corresponding argument of [pm.array.logspace](@ref logspace).
%>  \param[in]  logub   :   See the documentation of the corresponding argument of [pm.array.logspace](@ref logspace).
%>  \param[in]  logskip :   See the documentation of the corresponding argument of [pm.array.logspace](@ref logspace).
%>  \param[in]  base    :   See the documentation of the corresponding argument of [pm.array.logspace](@ref logspace).
%>
%>  \return
%>   `array`            :   The output vector of MATLAB real values containing
%>                          the set of logarithmically-spaced integer values in
%>                          the specified input range with the specified ``base``.
%>
%>  \interface{logspaceint}
%>  \code{.m}
%>
%>      array = pm.array.logspaceint(loglb, logub)
%>      array = pm.array.logspaceint(loglb, logub, logskip)
%>      array = pm.array.logspaceint(loglb, logub, [], base)
%>      array = pm.array.logspaceint(loglb, logub, logskip, base)
%>  \endcode
%>
%>  \example{logspaceint}
%>  \include{lineno} example/array/logspaceint/main.m
%>  \output{logspaceint}
%>  \include{lineno} example/array/logspaceint/main.out.m
%>
%>  \final{logspaceint}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 3:41 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function array = logspaceint(loglb, logub, logskip, base)
    if nargin < 3
        logskip = [];
    end
    if nargin < 4
        array = unique(round(pm.array.logspace(loglb, logub, logskip)));
    else
        array = unique(round(pm.array.logspace(loglb, logub, logskip, base)));
    end
end
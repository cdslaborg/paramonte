%>  \brief
%>  Return ``log(exp(logLarger) - exp(logSmaller))`` more accurately.
%>
%>  \warning
%>  The onus is on the user to ensure ``logSmaller < logLarger``.
%>
%>  \param[in]  logSmaller  :   The input scalar MATLAB ``real`` representing
%>                              the natural logarithm of the smaller value.
%>
%>  \param[in]  logLarger   :   The input scalar MATLAB ``real`` representing
%>                              the natural logarithm of the larger value.
%>
%>  \return
%>  ``logDiff``             :   The output scalar MATLAB ``real`` containing the
%>                              result of ``log(exp(logLarger) - exp(logSmaller))`` accurately.
%>
%>  \interface{logsubexp}
%>  \code{.m}
%>
%>      logDiff = pm.math.logsubexp(smaller, larger)
%>
%>  \endcode
%>
%>  \example{logsubexp}
%>  \include{lineno} example/math/logsubexp/main.m
%>  \output{logsubexp}
%>  \include{lineno} example/math/logsubexp/main.out.m
%>
%>  \final{logsubexp}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 9:45 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function logDiff = logsubexp(logSmaller, logLarger)
    logDiff = logLarger + log(-expm1(logSmaller - logLarger));
end
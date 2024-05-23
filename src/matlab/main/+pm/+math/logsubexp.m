%>  \brief
%>  Return ``log(exp(larger) - exp(smaller))`` more accurately.
%>
%>  \warning
%>  The onus is on the user to ensure ``smaller < larger``.
%>
%>  \param[in]  smaller :   The input scalar MATLAB real
%>                          representing the natural logarithm
%>                          of the smaller value.
%>
%>  \param[in]  smaller :   The input scalar MATLAB real
%>                          representing the natural logarithm
%>                          of the larger value.
%>
%>  \return
%>  val
%>      The output scalar MATLAB real containing the
%>      the result of ``log(exp(larger) - exp(smaller))`` accurately.
%>
%>  \interface{logsubexp}
%>  \code{.m}
%>
%>      val = pm.math.logsubexp(smaller, larger)
%>
%>  \endcode
%>  \example{logsubexp}
%>
%>      pm.math.logsubexp(log(400), log(1000))
%>
%>  \final{logsubexp}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 9:45 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function val = logsubexp(smaller, larger)
    val = larger + log(-expm1(smaller - larger));
end
%>  \brief
%>  Return a set of maximum ``sizemax`` unique integer spacings almost
%>  linearly spaced in the natural logarithmic scale between the
%>  specified ``start`` and ``stop`` of the output range.
%>
%>  \param[in]  start   :   The input scalar MATLAB whole-number (integer)
%>                          representing the starting point of the output range.
%>  \param[in]  stop    :   The input scalar MATLAB whole-number (integer)
%>                          representing the stopping point of the output range.
%>  \param[in]  sizemax :   The input scalar MATLAB whole-number (integer)
%>                          representing the maximum size of the output range.<br>
%>                          Due to rounding operation involved in creating the
%>                          output range, it is impossible to prespecify the
%>                          output range size, only the maximum.<br>
%>                          (**optional**, default = ``100``)
%>
%>  \return
%>  ``array``           :   The output vector of MATLAB real values containing
%>                          the set of naturally logarithmically-spaced integer
%>                          values in the specified input range.
%>
%>  \interface{logrange}
%>  \code{.m}
%>
%>      array = pm.array.logrange(start, stop)
%>      array = pm.array.logrange(start, stop, [])
%>      array = pm.array.logrange(start, stop, sizemax)
%>
%>  \endcode
%>
%>  \example{logrange}
%>  \include{lineno} example/array/logrange/main.m
%>  \output{logrange}
%>  \include{lineno} example/array/logrange/main.out.m
%>
%>  \final{logrange}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:29 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function array = logrange(start, stop, sizemax)
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
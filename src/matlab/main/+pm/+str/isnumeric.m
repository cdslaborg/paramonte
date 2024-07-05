%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the input string can be converted to a **scalar** number.<br>
%>
%>  \details
%>  The returned result is ``~isnan(str2double(str))``.<br>
%>  This is different from the result returned by
%>  the MATLAB intrinsic ``isnumeric()``.<br>
%>
%>  \param[in]  str :   The input scalar MATLAB string
%>                      whose conversion to numeric value is to be tested.<br>
%>
%>  \return
%>  ``itis``        :   The output scalar MATLAB logical that is ``true`` if and
%>                      only if the input ``str`` contains text that is convertible
%>                      to a **scalar** number, e.g., integer, real, complex.<br>
%>
%>  \interface{isnumeric}
%>  \code{.m}
%>
%>      itis = pm.str.isnumeric(str)
%>
%>  \endcode
%>
%>  \example{isnumeric}
%>  \include{lineno} example/str/isnumeric/main.m
%>  \output{isnumeric}
%>  \include{lineno} example/str/isnumeric/main.out.m
%>
%>  \final{isnumeric}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:36 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = isnumeric(str)
    itis = ~isnan(str2double(str));
end
%>  \brief
%>  Return the input scalar MATLAB string as doubly quoted string
%>  such that the first and last character of the output string are double-quotation
%>  marks while escaping any instances of **double quote** in the Fortran/MATLAB style,
%>  by duplicating every quotation mark within the string.<br>
%>
%>  \param[in]  str :   The input scalar MATLAB string to be doubly quoted.<br>
%>
%>  \return
%>  ``strQuoted``   :   The output scalar MATLAB string containing the doubly-quoted escaped input string.<br>
%>
%>  \interface{quote}
%>  \code{.m}
%>
%>      strQuoted = pm.str.quote(str)
%>
%>  \endcode
%>
%>  \example{quote}
%>  \include{lineno} example/str/quote/main.m
%>  \output{quote}
%>  \include{lineno} example/str/quote/main.out.m
%>
%>  \final{quote}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:40 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function strQuoted = quote(str)
    strQuoted = """" + strrep(str, """", """""") + """";
end

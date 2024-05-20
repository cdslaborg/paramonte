%
%   Return the input scalar MATLAB string as doubly quoted string
%   such that the first and last character of the output string are double-quotation
%   marks while escaping any instances of double quote in the Fortran/MATLAB style,
%   by duplicating every quotation mark within the string.
%
%       str
%
%           The input scalar MATLAB string to be doubly quoted.
%
%>  \return
%       strQuoted
%
%           The output scalar MATLAB string containing the doubly-quoted escaped input string.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       strQuoted = pm.str.quote(str)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function strQuoted = quote(str)
    strQuoted = """" + strrep(str, """", """""") + """";
end

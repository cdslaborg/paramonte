%>  \brief
%>  Return a vector of real values by splitting the
%>  input ``str`` string by the specified input ``sep``.<br>
%>
%>  \details
%>  This function is primarily used for splitting a
%>  string into a double vector of real numbers.<br>
%>
%>  \param[in]  str :   The input scalar MATLAB string containing from which the numbers must be retrieved
%>                      such as software version corresponding to major, minor, and patch versions.<br>
%>  \param[in]  sep :   The input scalar MATLAB string containing a string
%>                      to be used the separator in splitting.<br>
%>
%>  \return
%>  ``real``        :   The output MATLAB vector of size ``nsep + 1`` containing
%>                      the numbers retrieved from the input version string, where ``nsep``
%>                      represents the number of input ``sep`` found in the input string.<br>
%>                      If the specified input ``sep`` is empty ``[]``, the default value is used.<br>
%>                      (**optional**, default = ``"."``)
%>
%>  \interface{split2real}
%>  \code{.m}
%>
%>      real = pm.str.split2real(str)
%>      real = pm.str.split2real(str, [])
%>      real = pm.str.split2real(str, sep)
%>
%>  \endcode
%>
%>  \final{split2real}
%>
%>  \example{split2real}
%>  \include{lineno} example/str/split2real/main.m
%>  \output{split2real}
%>  \include{lineno} example/str/split2real/main.out.m
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:41 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function real = split2real(str, sep)
    if  nargin < 2
        sep = [];
    end
    if  isempty(sep)
        sep = ".";
    end
    real = str2double(string(strsplit(str, sep)));
end

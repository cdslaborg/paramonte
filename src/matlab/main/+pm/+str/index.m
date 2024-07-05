%>  \brief
%>  Return the starting index of the first occurrence
%>  of the input scalar MATLAB string ``pattern`` in
%>  the input scalar MATLAB string ``str``, otherwise,
%>  return ``0`` to indicate the lack of the ``pattern``
%>  in the input ``str``.<br>
%>
%>  \details
%>  This function partially replicates the functionality
%>  of the Fortran intrinsic function ``index()``.<br>
%>  This function uses the MATLAB intrinsic ``strfind()``
%>  to achieve the goal. However, unlike ``strfind()``
%>  it always returns a number such that the function
%>  can be directly used in string slicing.<br>
%>
%>  \param[in]  str     :   The input scalar MATLAB string to be searched
%>                          for the presence of the input pattern.<br>
%>  \param[in]  pattern :   The input scalar MATLAB string to be
%>                          searched for within the input ``str``.<br>
%>
%>  \return
%>  ``loc``             :   The output scalar MATLAB integer
%>                          containing the location of the first occurrence of the
%>                          input ``pattern`` in the input ``str`` or ``0`` if no such
%>                          pattern exists.<br>
%>
%>  \interface{index}
%>  \code{.m}
%>
%>      loc = pm.str.index(str, pattern)
%>
%>  \endcode
%>
%>  \example{index}
%>  \include{lineno} example/str/index/main.m
%>  \output{index}
%>  \include{lineno} example/str/index/main.out.m
%>
%>  \final{index}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:34 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function loc = index(str, pattern)
    loc = strfind(str, pattern);
    if isempty(loc)
        loc = 0;
    else
        loc = loc(1);
    end
end
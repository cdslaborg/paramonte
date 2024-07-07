%>  \brief
%>  Return a MATLAB string containing the ParaMonte MATLAB library authors.
%>
%>  \return
%>  ``str`` :   The output MATLAB string containing the ParaMonte MATLAB library authors.
%>
%>  \interface{authors}
%>  \code{.m}
%>
%>      str = pm.lib.authors();
%>
%>  \endcode
%>
%>  \example{authors}
%>  \include{lineno} example/lib/authors/main.m
%>  \output{authors}
%>  \include{lineno} example/lib/authors/main.out.m
%>
%>  \final{authors}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:59 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = authors()
    str = "The Computational Data Science Lab @ The University of Texas";
end
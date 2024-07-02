%>  \brief
%>  Return a  MATLAB string containing the ParaMonte MATLAB library acknowledgment statement.
%>
%>  \return
%>  A MATLAB string containing the ParaMonte MATLAB library acknowledgment statement.
%>
%>  \interface{credits}
%>  \code{.m}
%>
%>      result = pm.lib.credits()
%>
%>  \endcode
%>
%>  \example{credits}
%>  \include{lineno} example/lib/credits/main.m
%>  \output{credits}
%>  \include{lineno} example/lib/credits/main.out.m
%>
%>  \final{credits}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:10 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function result = credits()
    result = "Peter O'Donnell Fellowship / Texas Advanced Computing Center";
end
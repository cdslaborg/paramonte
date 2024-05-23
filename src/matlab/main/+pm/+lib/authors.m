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
%>  \final{authors}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:59 PM, University of Texas at Arlington<br>
function str = authors()
    str = "The Computational Data Science Lab @ The University of Texas";
end
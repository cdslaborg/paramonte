%>  \brief
%>  Return a scalar MATLAB string containing ``4``
%>  blank characters equivalent to a tab character.
%>
%>  \details
%>  This function primarily exists to bring
%>  consistency to the ParaMonte library IO tasks.
%>
%>  \return
%>  `str`   :   The output scalar MATLAB string containing ``4``
%>              blank characters equivalent to a tab character.
%>
%>  \interface{tab}
%>  \code{.m}
%>
%>      str = pm.io.tab()
%>
%>  \endcode
%>
%>  \final{tab}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:58 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = tab()
    str = "    ";
end
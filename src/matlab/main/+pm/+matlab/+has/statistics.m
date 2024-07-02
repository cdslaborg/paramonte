%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB Statistics Toolbox.
%>
%>  \details
%>  This function searches the MATLAB license
%>  for an installation of the Toolbox.
%>
%>  \return
%>  ``hasit``   :   The output scalar MATLAB logical that is ``true`` if and
%>                  only if the current installation of MATLAB contains
%>                  the required MATLAB Toolbox.
%>
%>  \interface{statistics}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.statistics();
%>
%>  \endcode
%>
%>  \example{statistics}
%>  \include{lineno} example/matlab/has/main.m
%>  \output{statistics}
%>  \include{lineno} example/matlab/has/main.out.m
%>
%>  \final{statistics}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:33 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = statistics()
    hasit = license('test', 'Statistics_Toolbox');
end
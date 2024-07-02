%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB mapping Toolbox.
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
%>  \interface{map}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.map();
%>
%>  \endcode
%>
%>  \example{map}
%>  \include{lineno} example/matlab/has/main.m
%>  \output{map}
%>  \include{lineno} example/matlab/has/main.out.m
%>
%>  \final{map}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:33 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = map()
    hasit = license('test', 'MAP_Toolbox');
end
%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB Econometrics Toolbox.
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
%>  \interface{econometrics}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.econometrics();
%>
%>  \endcode
%>
%>  \example{econometrics}
%>  \include{lineno} example/matlab/has/main.m
%>  \output{econometrics}
%>  \include{lineno} example/matlab/has/main.out.m
%>
%>  \final{econometrics}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:00 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = econometrics()
    hasit = license('test', 'Econometrics_Toolbox');
end
%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB Medical Imaging Toolbox.
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
%>  \interface{medimage}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.medimage();
%>
%>  \endcode
%>
%>  \example{medimage}
%>  \include{lineno} example/matlab/has/main.m
%>  \output{medimage}
%>  \include{lineno} example/matlab/has/main.out.m
%>
%>  \final{medimage}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:33 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = medimage()
    hasit = license('test', 'Medical_Imaging_Toolbox');
end
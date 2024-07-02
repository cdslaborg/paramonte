%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB Bioinformatics Toolbox.
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
%>  \interface{bioinformatics}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.bioinformatics();
%>
%>  \endcode
%>
%>  \example{bioinformatics}
%>  \include{lineno} example/matlab/has/main.m
%>  \output{bioinformatics}
%>  \include{lineno} example/matlab/has/main.out.m
%>
%>  \final{bioinformatics}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 9:49 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = bioinformatics()
    hasit = license('test', 'Bioinformatics_Toolbox');
end
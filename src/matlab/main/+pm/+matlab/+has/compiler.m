%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB Compiler Toolbox.
%>
%>  \details
%>  This function searches the MATLAB license
%>  for an installation of the Toolbox.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  hasit
%>      The output scalar MATLAB logical that is ``true`` if and
%>      only if the current installation of MATLAB contains
%>      the required MATLAB Toolbox.
%>
%>  \interface{compiler}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.compiler();
%>
%>  \endcode
%>  \final{compiler}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 9:51 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = compiler()
    hasit = license('test', 'Compiler');
end
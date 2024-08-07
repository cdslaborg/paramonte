%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the current installation of MATLAB contains
%>  the MATLAB Parallel Computing Toolbox.
%>
%>  \details
%>  This function searches the MATLAB license
%>  for an installation of the Parallel Computing Toolbox.
%>  If the search fails, a parallel code section will be tested.
%>
%>  \return
%>  ``hasit``   :   The output scalar MATLAB logical that is ``true`` if and
%>                  only if the current installation of MATLAB contains
%>                  the required MATLAB Toolbox.
%>
%>  \interface{parallel}
%>  \code{.m}
%>
%>      hasit = pm.matlab.has.parallel();
%>
%>  \endcode
%>
%>  \example{parallel}
%>  \include{lineno} example/matlab/has/main.m
%>  \output{parallel}
%>  \include{lineno} example/matlab/has/main.out.m
%>
%>  \final{parallel}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:33 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hasit = parallel()
    hasit = license('test', 'Distrib_Computing_Toolbox');
    if ~hasit
        try
            delete(gcp('nocreate'));
            hasit = true;
        end
    end
end
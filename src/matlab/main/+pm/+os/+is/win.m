%>  \brief
%>  Return ``true`` if the current OS is Windows.
%>
%>  \return
%>  ``itis``    :   The output MATLAB logical scalar value that is ``true``
%>                  if and only if the OS is Windows, otherwise ``false``.
%>
%>  \interface{win}
%>  \code{.m}
%>
%>      itis = pm.os.is.win()
%>
%>  \endcode
%>
%>  \example{win}
%>  \include{lineno} example/os/is/main.m
%>  \output{win}
%>  \include{lineno} example/os/is/main.out.m
%>
%>  \final{win}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:47 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = win()
    itis = ispc;
end
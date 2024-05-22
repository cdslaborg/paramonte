%>  \brief
%>  Return ``true`` if the current OS is Linux.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `itis`  :   The output MATLAB logical scalar value that is ``true``
%>              if and only if the OS is Linux, otherwise ``false``.
%>
%>  \interface{lin}
%>  \code{.m}
%>
%>      itis = pm.os.is.lin()
%>
%>  \endcode
%>  \final{lin}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:45 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = lin()
    itis = isunix && ~ismac;
end
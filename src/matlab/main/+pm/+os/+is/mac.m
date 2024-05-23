%>  \brief
%>  Return ``true`` if the current OS is MacOS (Darwin).
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `itis`  :   The output MATLAB logical scalar value that is ``true`` 
%>              if and only if the OS is MacOS (Darwin), otherwise ``false``.
%>
%>  \interface{mac}
%>  \code{.m}
%>      itis = pm.os.is.mac()
%>
%>  \endcode
%>
%>  \final{mac}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:46 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function result = mac()
    result = ismac;
end
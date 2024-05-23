%>  \brief
%>  Return a list of MATLAB strings containing the names of
%>  all supported compiler suites (vendor names) used for
%>  building the ParaMonte MATLAB shared libraries.
%>
%>  \devnote
%>  The CS names within this function must be
%>  regularly updated with the latest CS names
%>  available in the ParaMonte installation guide.
%>
%>  \params `None`
%>
%>  \return
%>  `typelist`  :   A MATLAB string list containing:
%>                  ``["intel", "gnu"]``
%>
%>  \interface{clstypes}
%>  \code{.m}
%>
%>      typelist = pm.lib.clstypes()
%>
%>  \endcode
%>  \final{clstypes}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:07 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function typelist = clstypes()
    typelist = ["intel", "gnu"];
end
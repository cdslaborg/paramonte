%>  \brief
%>  Return a list of MATLAB strings containing the names of
%>  all currently possible builds of the ParaMonte MATLAB
%>  shared libraries.
%>
%>  \return
%>  ``typelist``    :   The output MATLAB string list containing the 
%>                      value ``["native", "tuned", "ipo", "release", "testing", "debug"]``.
%>
%>  \interface{bldtypes}
%>  \code{.m}
%>
%>      typelist = pm.lib.bldtypes();
%>
%>  \endcode
%>
%>  \devnote
%>  The build names within this function must be
%>  regularly updated with the latest build names
%>  available in the ParaMonte installation guide.
%>
%>  \final{bldtypes}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:01 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function typelist = bldtypes()
    typelist = ["native", "tuned", "ipo", "release", "testing", "debug"];
end
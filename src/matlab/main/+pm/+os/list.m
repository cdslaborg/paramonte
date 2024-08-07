%>  \brief
%>  Return a list of MATLAB strings containing the names of
%>  OS platforms supported by the ParaMonte MATLAB library.
%>
%>  \return
%>  ``names``   :   The output MATLAB string list containing:
%>                  ``["Windows", "Linux", "Darwin"]``.
%>
%>  \interface{list}
%>  \code{.m}
%>
%>      names = pm.os.list()
%>
%>  \endcode
%>
%>  \example{list}
%>  \include{lineno} example/os/list/main.m
%>  \output{list}
%>  \include{lineno} example/os/list/main.out.m
%>
%>  \final{list}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:49 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function names = list()
    names = ["Windows", "Linux", "Darwin"];
end
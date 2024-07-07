%>  \brief
%>  Return a list of MATLAB strings containing the lower-case names
%>  of all OS platforms supported by the ParaMonte MATLAB library.
%>
%>  \return
%>  ``names``   :   The output MATLAB string list containing:
%>                  ``["windows", "linux", "darwin"]``.
%>
%>  \interface{listl}
%>  \code{.m}
%>
%>      names = pm.os.listl()
%>
%>  \endcode
%>
%>  \example{listl}
%>  \include{lineno} example/os/listl/main.m
%>  \output{listl}
%>  \include{lineno} example/os/listl/main.out.m
%>
%>  \final{listl}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:51 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function names = listl()
    names = lower(pm.os.list());
end
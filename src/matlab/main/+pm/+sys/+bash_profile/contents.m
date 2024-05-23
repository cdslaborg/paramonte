%>  \brief
%>  Return the contents of the ``.bash_profile`` file 
%>  in the system home folder as a scalar MATLAB string.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `str`   :   The output scalar MATLAB string containing the contents
%>              of the ``.bash_profile`` file if it exists or an empty string.
%>
%>  \interface{contents}
%>  \code{.m}
%>
%>      str = pm.sys.bash_profile.contents()
%>
%>  \endcode
%>  \final{contents}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:57 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = contents()
    path = fullfile(pm.sys.path.home(), ".bash_profile");
    if  isfile(path)
        str = string(fileread(path));
    else
        str = "";
    end
end
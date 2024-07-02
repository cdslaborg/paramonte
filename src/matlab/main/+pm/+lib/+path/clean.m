%>  \brief
%>  Remove all paths that contain the ParaMonte
%>  `lib` directory from the MATLAB ``path`` variable.
%>
%>  \details
%>  This action is particularly vital for correct handling
%>  of MATLAB MEX files for different build configurations.
%>
%>  \interface{clean}
%>  \code{.m}
%>
%>      pm.lib.path.clean()
%>
%>  \endcode
%>
%>  \final{clean}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:31 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function clean()
    libpath = pm.lib.path.lib();
    pathlist = strsplit(string(path), pathsep);
    for ipath = 1 : length(pathlist)
        if  contains(pathlist{ipath}, libpath)
            rmpath(pathlist{ipath});
        end
    end
end
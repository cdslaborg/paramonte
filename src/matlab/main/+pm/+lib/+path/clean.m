%
%   Remove all paths that contain the ParaMonte
%   `lib` directory from the MATLAB ``path`` variable.
%
%   This action is particularly vital for correct handling
%   of MATLAB MEX files for different build configurations.
%
%       None
%
%>  \return
%       None
%
%   Interface
%   ---------
%
%       pm.lib.path.clean()
%
%>  \final{}
%>
%>  \author
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
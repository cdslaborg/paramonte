%>  \brief
%>  Return a MATLAB string containing the
%>  absolute path to the system home directory.<br>
%>
%>  \return
%>  ``path``    :   A MATLAB string containing the
%>                  absolute path to the system home directory.<br>
%>
%>  \interface{home}
%>  \code{.m}
%>
%>      path = pm.sys.path.home()
%>
%>  \endcode
%>
%>  \example{home}
%>  \include{lineno} example/sys/path/home/main.m
%>  \output{home}
%>  \include{lineno} example/sys/path/home/main.out.m
%>
%>  \final{home}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:25 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function path = home()
    persistent homePathPersistent
    %freshRequested = false;
    %if nargin==0
    %    if isempty(homePathPersistent); freshRequested = true; end
    %elseif nargin==1 && ( strcmpi(string(varargin{1}),"fresh") || strcmpi(string(varargin{1}),"new") || strcmpi(string(varargin{1}),"reset") )
    %    freshRequested = true;
    %else
    %    error( "home() takes at most one argument of the following values: ""new"", ""fresh"", ""reset"", all with the same meaning." );
    %end
    %if freshRequested
    if ispc
        [failed, homePathPersistent] = system("echo %HOMEPATH%");
        if failed
            error   ( newline ...
                    + "Failed to capture the path to the home directory of the Windows system." + newline ...
                    + "This is highlt unusual and likely a low-level MATLAB problem." + newline ...
                    + newline ...
                    );
        else
            homePathPersistent = strtrim(homePathPersistent);
        end
    else
        homePathPersistent = strtrim(pm.sys.path.abs("~"));
    end
    %end
    path = string(homePathPersistent);
end
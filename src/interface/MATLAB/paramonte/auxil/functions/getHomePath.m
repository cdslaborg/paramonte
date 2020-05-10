function homePath = getHomePath(varargin)
    freshRequested = false;
    persistent homePathPersistent
    if nargin==0
        if isempty(homePathPersistent); freshRequested = true; end
    elseif nargin==1 && ( strcmpi(string(varargin{1}),"fresh") || strcmpi(string(varargin{1}),"new") || strcmpi(string(varargin{1}),"reset") )
        freshRequested = true;
    else
        error( "getHomePath() takes at most one argument of the following values: ""new"", ""fresh"", ""reset"", all with the same meaning." );
    end
    if freshRequested
        if ispc
            [errStat,homePathPersistent] = system("echo %HOMEPATH%");
            if errStat  
                error   ( newline ...
                        + "Failed to capture the path to the home directory of the Windows system. This is highlt unusual and likely a low-level MATLAB problem." ...
                        + newline ...
                        );
            else
                homePathPersistent = strtrim(homePathPersistent);
            end
        else
            homePathPersistent = strtrim(getFullPath("~"));
        end
    end
    homePath = string(homePathPersistent);
end
function mexpaths = mex(mexname, config)
    %
    %   Return the vector of MATLAB strings containing the
    %   path(s) to the ParaMonte library mex shared files.
    %
    %   Parameters
    %   ----------
    %
    %       mexname
    %
    %           The input scalar MATLAB string containing
    %           the (partial) name of a MEX file name without MEX extension.
    %           (**optional**. If missing, all detected paths will be returned.)
    %
    %       config
    %
    %           The input scalar (or vector of) MATLAB string(s) containing the
    %           (set of) keyword(s) to match in the identified MEX files paths.
    %           (**optional**. If missing, all detected paths will be returned.)
    %
    %   Returns
    %   -------
    %
    %       mexpaths
    %
    %           The output vector of MATLAB strings containing the
    %           path(s) to the ParaMonte library MEX shared files.
    %           The output ``mexpaths`` will be an empty list
    %           if no such paths are found.
    %
    %   Interface
    %   ---------
    %
    %       mexpaths = pm.lib.path.mex()
    %       mexpaths = pm.lib.path.mex([])
    %       mexpaths = pm.lib.path.mex([], [])
    %       mexpaths = pm.lib.path.mex(mexname)
    %       mexpaths = pm.lib.path.mex([], config)
    %       mexpaths = pm.lib.path.mex(mexname, [])
    %       mexpaths = pm.lib.path.mex(mexname, config)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if  nargin < 2
        config = [];
    end
    if  nargin < 1
        mexname = [];
    end
    if ~isempty(config)
        pm.introspection.verify(config, "string", 255, "config");
        %configList = strsplit(config, ":");
        configList = string(config);
    else
        configList = [];
    end
    if ~isempty(mexname)
        pm.introspection.verify(mexname, "string", 1, "mexname");
    else
        mexname = "";
    end
    libdir = pm.lib.path.lib();
    mexpaths = pm.sys.path.glob(fullfile(libdir, "**" + mexname + "*" + mexext));
    if ~isempty(configList) && ~isempty(mexpaths)
        libdirlen = length(libdir{1});
        indices = [];
        for i = 1 : length(mexpaths)
            path = mexpaths{i};
            matching = true;
            for j = 1 : length(configList)
                key = configList{j};
                if ~contains(path(libdirlen + 1 : end), key)
                    matching = false;
                    break;
                end
            end
            if matching
                indices = [indices, i];
            end
        end
        mexpaths = mexpaths(indices);
    end
end
function clean()
    %
    %   Remove all paths that contain the ParaMonte
    %   `lib` directory from the MATLAB ``path`` variable.
    %
    %   This action is particularly vital for correct handling
    %   of MATLAB MEX files for different build configurations.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       None
    %
    %   Interface
    %   ---------
    %
    %       pm.lib.path.clean()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    libpath = pm.lib.path.lib();
    pathlist = strsplit(string(path), pathsep);
    for ipath = 1 : length(pathlist)
        if  contains(pathlist{ipath}, libpath)
            rmpath(pathlist{ipath});
        end
    end
end
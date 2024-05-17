function path = lib(key)
    %
    %   Return a scalar MATLAB string containing the path to
    %   the ``lib`` directory of the ParaMonte library package.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       path
    %
    %           The output scalar MATLAB string containing the path to
    %           the ``lib`` directory of the ParaMonte library package.
    %
    %   Interface
    %   ---------
    %
    %       path = pm.lib.path.lib()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    persistent path_persistent
    if isempty(path_persistent)
        path_persistent = string(fullfile(pm.lib.path.root(), "lib"));
    end
    path = path_persistent;
end
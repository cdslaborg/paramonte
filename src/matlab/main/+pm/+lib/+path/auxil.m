function path = auxil()
    %
    %   Return a scalar MATLAB string containing the path to
    %   the ``auxil`` directory of the ParaMonte library package.
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
    %           the ``auxil`` directory of the ParaMonte library package.
    %
    %   Interface
    %   ---------
    %
    %       path = pm.lib.path.auxil()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    persistent path_persistent
    if isempty(path_persistent)
        path_persistent = string(fullfile(pm.lib.path.root(), "auxil"));
    end
    path = path_persistent;
end
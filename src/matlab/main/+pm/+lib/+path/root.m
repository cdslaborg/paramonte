function path = root()
    %
    %   Return a scalar MATLAB string containing the
    %   root directory of the ParaMonte library package.
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
    %           A scalar MATLAB string containing the root
    %           directory of the ParaMonte library package.
    %
    %   Interface
    %   ---------
    %
    %       path = pm.lib.path.root()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    persistent path_persistent
    if isempty(path_persistent)
        [path_persistent, ~, ~] = fileparts(mfilename('fullpath'));
        path_persistent = string(pm.sys.path.abs(fullfile(path_persistent, "..", "..")));
    end
    path = path_persistent;
    %addpath(genpath(tree.root),'-begin');
end
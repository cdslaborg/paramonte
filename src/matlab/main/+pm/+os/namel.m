function result = namel()
    %
    %   Return a MATLAB string containing the lower-case name of the current OS.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       str
    %
    %           The output MATLAB string containing either:
    %
    %               ``"linux"`` if the OS is Linux.
    %               ``"windows"`` if the OS is Windows.
    %               ``"darwin"`` if the OS is macOS (Darwin).
    %
    %   Interface
    %   ---------
    %
    %       str = pm.os.namel()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if ispc
        result = "windows";
    elseif ismac
        result = "darwin";
    elseif isunix
        result = "linux";
    end
end
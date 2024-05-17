function str = name()
    %
    %   Return a MATLAB string containing the name of the current OS.
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
    %           The output MATLAB string that is either:
    %
    %               ``"Linux"`` if the OS is Linux.
    %               ``"Windows"`` if the OS is Windows.
    %               ``"Darwin"`` if the OS is macOS (Darwin).
    %
    %   Interface
    %   ---------
    %
    %       str = pm.os.name()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if ispc
        str = "Windows";
    elseif ismac
        str = "Darwin";
    elseif isunix
        str = "Linux";
    end
end
function str = tab()
    %
    %   Return a scalar MATLAB string containing ``4``
    %   blank characters equivalent to a tab character.
    %
    %   This function primarily exists to bring
    %   consistency to the ParaMonte library IO tasks.
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
    %           The output scalar MATLAB string containing ``4``
    %           blank characters equivalent to a tab character.
    %
    %   Interface
    %   ---------
    %
    %       pm.io.tab()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    str = "    ";
end
function itis = win()
    %
    %   Return ``true`` if the current OS is Windows.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       itis
    %
    %           The output MATLAB logical scalar value that is ``true``
    %           if and only if the OS is Windows, otherwise ``false``.
    %
    %   Interface
    %   ---------
    %
    %       itis = pm.os.is.win()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    itis = ispc;
end
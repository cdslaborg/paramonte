function itis = lin()
    %
    %   Return ``true`` if the current OS is Linux.
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
    %           if and only if the OS is Linux, otherwise ``false``.
    %
    %   Interface
    %   ---------
    %
    %       itis = pm.os.is.lin()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    itis = isunix && ~ismac;
end
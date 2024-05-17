function result = mac()
    %
    %   Return ``true`` if the current OS is MacOS (Darwin).
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
    %           if and only if the OS is MacOS (Darwin), otherwise ``false``.
    %
    %   Interface
    %   ---------
    %
    %       itis = pm.os.is.mac()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    result = ismac;
end
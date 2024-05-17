function itis = iscmd()
    %
    %   Return a scalar MATLAB logical ``true`` if and only if the
    %   the MATLAB binary is being called from the shell command-line,
    %   in ``batch`` mode (without Graphical User Interface (GUI) interface).
    %   Otherwise, return ``false`` implying that MATLAB GUI engine is active.
    %
    %   This functionality is important for certain library
    %   features that require MATLAB to be called from the
    %   command-line, e.g., MPI-parallel functionalities.
    %
    %   \warning
    %
    %       This function relies on functionalities
    %       that are supported only in MATLAB > 2019a.
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
    %           The output scalar MATLAB logical that is ``true``
    %           if and only if ParaMonte MATLAB library is being called
    %           from MATLAB Graphical User Interface (GUI), otherwise, ``false``.
    %
    %   Interface
    %   ---------
    %
    %       itis = pm.matlab.iscmd()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    try
        %   This error catching can be removed once
        %   no MATLAB < R2019a is expected to call it.
        itis = batchStartupOptionUsed; % batchStartupOptionUsed is introduced in R2019a and not supported in older versions of MATLAB
    catch
        try
            itis = ~pm.matlab.isgui();
        catch
            itis = isdeployed();
        end
    end
end
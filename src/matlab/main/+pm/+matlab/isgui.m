function itis = isgui()
    %
    %   Return a scalar MATLAB logical ``true`` if and only if the
    %   current environment is the MATLAB Graphical User Interface (GUI),
    %   otherwise, return ``false`` implying that MATLAB engine is active
    %   without the GUI (e.g., called from the system shell command line).
    %
    %   This functionality is important for certain
    %   library features that require MATLAB GUI.
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
    %       itis = pm.matlab.isgui()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    persistent itis_persistent
    if isempty(itis_persistent)
        try
            % usejava() is sensitive to char vs. string input arguments. always input char.
            itis_persistent = usejava('desktop') && feature('ShowFigureWindows') && ~batchStartupOptionUsed; % batchStartupOptionUsed is introduced in R2019a and not supported in older versions of MATLAB
        catch
            itis_persistent = usejava('desktop') && feature('ShowFigureWindows');
        end
    end
    itis = itis_persistent;
end
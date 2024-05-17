function typelist = clstypes()
    %
    %   Return a list of MATLAB strings containing the names of
    %   all supported compiler suites (vendor names) used for
    %   building the ParaMonte MATLAB shared libraries.
    %
    %   \devnote
    %
    %       The CS names within this function must be
    %       regularly updated with the latest CS names
    %       available in the ParaMonte installation guide.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       typelist
    %
    %           A MATLAB string list containing:
    %
    %               ``["intel", "gnu"]``
    %
    %   Interface
    %   ---------
    %
    %       typelist = pm.lib.clstypes()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    typelist = ["intel", "gnu"];
end
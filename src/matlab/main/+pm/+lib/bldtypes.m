function typelist = bldtypes()
    %
    %   Return a list of MATLAB strings containing the names of
    %   all currently possible builds of the ParaMonte MATLAB
    %   shared libraries.
    %
    %   \devnote
    %
    %       The build names within this function must be
    %       regularly updated with the latest build names
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
    %               ``["native", "tuned", "ipo", "release", "testing", "debug"]``
    %
    %   Interface
    %   ---------
    %
    %       typelist = pm.lib.bldtypes()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    typelist = ["native", "tuned", "ipo", "release", "testing", "debug"];
end
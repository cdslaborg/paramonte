function names = listl()
    %
    %   Return a list of MATLAB strings containing the lower-case names
    %   of all OS platforms supported by the ParaMonte MATLAB library.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       names
    %
    %           The output MATLAB string list containing:
    %
    %           ``["windows", "linux", "darwin"]``
    %
    %   Interface
    %   ---------
    %
    %       names = pm.os.listl()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    names = lower(pm.os.list());
end
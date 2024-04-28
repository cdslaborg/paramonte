function names = list()
    %
    %   Return a list of MATLAB strings containing the names of
    %   OS platforms supported by the ParaMonte MATLAB library.
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
    %               ``["Windows", "Linux", "Darwin"]``
    %
    %   Interface
    %   ---------
    %
    %       names = pm.os.list()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    names = ["Windows", "Linux", "Darwin"];
end
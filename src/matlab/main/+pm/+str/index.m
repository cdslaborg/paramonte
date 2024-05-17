function loc = index(str, pattern)
    %
    %   Return the starting index of the first occurrence
    %   of the input scalar MATLAB string ``pattern`` in
    %   the input scalar MATLAB string ``str``, otherwise,
    %   return ``0`` to indicate the lack of the ``pattern``
    %   in the input ``str``.
    %
    %   This function partially replicates the functionality
    %   of the Fortran intrinsic function ``index()``.
    %
    %   This function uses the MATLAB intrinsic ``strfind()``
    %   to achieve the goal. However, unlike ``strfind()``
    %   it always returns a number such that the function
    %   can be directly used in string slicing.
    %
    %   Parameters
    %   ----------
    %
    %       str
    %
    %           The input scalar MATLAB string to be searched
    %           for the presence of the input pattern.
    %
    %       pattern
    %
    %           The input scalar MATLAB string to be
    %           searched for within the input ``str``.
    %
    %   Returns
    %   -------
    %
    %       loc
    %
    %           The output scalar MATLAB integer
    %           containing the location of the first occurrence of the
    %           input ``pattern`` in the input ``str`` or ``0`` if no such
    %           pattern exists.
    %
    %   Interface
    %   ---------
    %
    %       loc = pm.str.index(str, pattern)
    %
    %   Example
    %   -------
    %
    %       loc = pm.str.index("paramonte", "") % 1
    %       loc = pm.str.index("paramonte", "M") % 0
    %       loc = pm.str.index("paramonte", "mont") % 5
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    loc = strfind(str, pattern);
    if isempty(loc)
        loc = 0;
    else
        loc = loc(1);
    end
end
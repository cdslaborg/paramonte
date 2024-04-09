function indx = index(array)
    %
    %   Given ``array(1 : n)``, return the array index ``indx(1:n)``
    %   such that ``array(indx)`` is in ascending order.
    %
    %   This is a simple convenience wrapper around
    %   the MATLAB intrinsic function ``sort()``.
    %
    %   \note
    %
    %       The output ``indx`` is NOT the rank
    %       of the elements of the input array.S
    %
    %   Parameters
    %   ----------
    %
    %       array
    %
    %           The input MATLAB vector of sortable values
    %           (that can be passed to the MATLAB intrinsic function ``sort()``).
    %
    %   Returns
    %   -------
    %
    %       indx
    %
    %           The output MATLAB integer vector of the same size as the
    %           input ``array`` such that ``array(indx)`` is in ascending order.
    %
    %   Interface
    %   ---------
    %
    %       indx = pm.sort.index(array)
    %
    %   Example
    %   -------
    %
    %       indx = pm.sort.index([2.5, 4, 1.5, -3])
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    [~, indx] = sort(array);
end
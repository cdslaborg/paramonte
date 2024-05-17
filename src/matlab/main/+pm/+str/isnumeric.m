function itis = isnumeric(str)
    %
    %   Return a scalar MATLAB logical that is ``true`` if and
    %   only if the input string can be converted to a number.
    %
    %   The returned result is ``~isnan(str2double(str))``.
    %   This is different from the result returned by
    %   the MATLAB intrinsic ``isnumeric()``.
    %
    %   Parameters
    %   ----------
    %
    %       str
    %
    %           The input scalar MATLAB string
    %           whose conversion to numeric value is to be tested.
    %
    %   Returns
    %   -------
    %
    %       itis
    %
    %           The output scalar MATLAB logical that is ``true`` if and
    %           only if the input ``str`` contains text that is convertible
    %           to number(s), e.g., integer, real, complex.
    %
    %   Interface
    %   ---------
    %
    %       itis = pm.str.isnumeric(str)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    itis = ~isnan(str2double(str));
end
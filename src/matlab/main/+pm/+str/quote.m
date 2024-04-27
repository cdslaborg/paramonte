function strQuoted = quote(str)
    %
    %   Return the input scalar MATLAB string as doubly quoted string
    %   such that the first and last character of the output string are double-quotation
    %   marks while escaping any instances of double quote in the Fortran/MATLAB style,
    %   by duplicating every quotation mark within the string.
    %
    %   Parameters
    %   ----------
    %
    %       str
    %
    %           The input scalar MATLAB string to be doubly quoted.
    %
    %   Returns
    %   -------
    %
    %       strQuoted
    %
    %           The output scalar MATLAB string containing the doubly-quoted escaped input string.
    %
    %   Interface
    %   ---------
    %
    %       strQuoted = pm.str.quote(str)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    strQuoted = """" + strrep(str, """", """""") + """";
end

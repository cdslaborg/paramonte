function real = split2real(str)
    %
    %   Return an integer triplet list from the input dot-separated string.
    %   This function is primarily used for splitting a dot-separated version string
    %   into a double vector of whole-number version.
    %
    %   Parameters
    %   ----------
    %
    %       str
    %
    %           The input MATLAB string containing a dot-separated string,
    %           such as software version corresponding to major, minor, and patch versions.
    %
    %   Returns
    %   -------
    %
    %       real
    %
    %           The output MATLAB vector of size ``ndot + 1`` containing the
    %           whole-number versions retrieved from the input version string,
    %           where ``ndot`` represents the number of dots in the input string.
    %
    %   Interface
    %   ---------
    %
    %       real = pm.str.split2real()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    real = str2double(string(strsplit(str, ".")));
end

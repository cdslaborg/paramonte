%
%   Return an integer triplet list from the input dot-separated string.
%   This function is primarily used for splitting a dot-separated version string
%   into a double vector of whole-number version.
%
%       str
%
%           The input MATLAB string containing a dot-separated string,
%           such as software version corresponding to major, minor, and patch versions.
%
%>  \return
%       real
%
%           The output MATLAB vector of size ``ndot + 1`` containing the
%           whole-number versions retrieved from the input version string,
%           where ``ndot`` represents the number of dots in the input string.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       real = pm.str.split2real()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function real = split2real(str)
    real = str2double(string(strsplit(str, ".")));
end

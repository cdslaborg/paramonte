%
%   Return the integers ``ncol`` and ``nrow`` such that their sum
%   (half the circumference of the resulting rectangle) is minimum while
%   their product is larger than or equal the input integer-valued ``area``.
%
%   This function is particularly useful for computing
%   the number of axes columns and rows in a single figure
%   that yield the most square-like rectangle.
%
%       area
%
%           The input scalar MATLAB whole number
%           representing the number of axes (subplots)
%           to arrange in a single figure.
%
%>  \return
%       nrow
%
%           The output scalar MATLAB whole-number
%           representing the number of rows of the rectangle.
%
%       ncol
%
%           The output scalar MATLAB whole-number
%           representing the number of columns of the rectangle.
%
%   Interface
%   ---------
%
%       [nrow, ncol] = pm.math.mincirc(area); 
%
%   Example
%   -------
%
%       [nrow, ncol] = pm.math.mincirc(5); % 3, 2
%       [nrow, ncol] = pm.math.mincirc(21); % 7, 3
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [nrow, ncol] = mincirc(area)
    ncol = 1:ceil(sqrt(area));
    nrow = ceil(area ./ ncol);
    [~, minloc] = min(nrow + ncol);
    nrow = nrow(minloc);
    ncol = ncol(minloc);
end
%
%   Return the input hashmap with the new (key, val) appended.
%   only if the input pair does not exist in the input hashmap.
%   If the pair exists, replace it with the new value.
%
%       key
%
%           The input scalar MATLAB string
%           containing the key to add to the input pair list.
%
%       val
%
%           The input value to be appended to the 
%           input hashmap after the input ``key``.
%
%       hashmap
%
%           The input cell array of even number of elements
%           containing the ``(key, val)`` pairs of the hashmap
%           in sequence as element of the cell array.
%
%>  \return
%       hashnew
%
%           The output cell array of even number of elements
%           containing the input pair list and ``(key, val)`` pair.
%           If the input ``key`` exists in the input ``hashmap``,
%           the input ``(key, val)`` pair will not be added.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       hashnew = pm.matlab.hashmap.repKeyVal(key, val, hashmap)
%
%   Example
%   -------
%
%       hashmap = {"key1", 1, "key2", "val2"};
%       hashnew = pm.matlab.hashmap.repKeyVal("key2", 2, hashmap)
%       hashnew = pm.matlab.hashmap.repKeyVal("key3", true, hashmap)
%       hashnew = pm.matlab.hashmap.repKeyVal("key3", false, hashmap)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hashnew = repKeyVal(key, val, hashmap)
    hashnew = hashmap;
    vararginLen = length(hashnew);
    keyString = string(key);
    failed = true;
    for i = 1 : vararginLen - 1
        try
            if strcmpi(string(hashnew{i}), keyString)
                hashnew{i + 1} = val;
                failed = false;
                break;
            end
        catch
            continue
        end
    end
    if  failed
        hashnew = {hashnew{:}, key, val};
    end
end
function hashnew = addKeyVal(key, val, hashmap)
    %
    %   Return the input hashmap with the new (key, val) appended,
    %   only if the input pair does not exist in the input hashmap.
    %
    %   Parameters
    %   ----------
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
    %           containing a set of ``(key, val)`` pairs of the hashmap.
    %
    %   Returns
    %   -------
    %
    %       hashnew
    %
    %           The output cell array of even number of elements
    %           containing the input pair list and ``(key, val)`` pair.
    %           If the input ``key`` exists in the input ``hashmap``,
    %           the input ``(key, val)`` pair will not be added.
    %
    %   Interface
    %   ---------
    %
    %       hashnew = pm.matlab.hashmap.addKeyVal(key, val, hashmap)
    %
    %   Example
    %   -------
    %
    %       hashmap = {"key1", 1, "key2", "val2"};
    %       hashnew = pm.matlab.hashmap.addKeyVal("key3", false, {})
    %       hashnew = pm.matlab.hashmap.addKeyVal("key2", "val2", hashmap) % = hashmap
    %       hashnew = pm.matlab.hashmap.addKeyVal("key3", false, hashmap)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    hashnew = hashmap;
    if  2 < nargin && 0 < length(hashnew)
        [currentVal, failed] = pm.matlab.hashmap.getVal(key, hashmap);
        if  failed
            hashnew = {hashnew{:}, key, val};
        end
    else
        hashnew = {key, val}; % XXX is this correct behavior?
    end
end
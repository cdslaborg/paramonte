function [val, failed] = getVal(key, hashmap)
    %
    %   Return the value corresponding to the input ``key``
    %   in the input ``hashmap`` cell array.
    %
    %   Parameters
    %   ----------
    %
    %       key
    %
    %           The input scalar MATLAB string
    %           containing the key to search for in the input pair list.
    %
    %       hashmap
    %
    %           The input cell array of even number of elements
    %           containing the ``(key, val)`` pairs of the hashmap
    %           in sequence as element of the cell array.
    %
    %   Returns
    %   -------
    %
    %       val
    %
    %           The output MATLAB object stored in the element of the
    %           input cell array, whose ``key`` is given as the input.
    %
    %       failed
    %
    %           The output scalar MATLAB logical that is ``true``
    %           if and only if the input ``key`` exists in the input ``hashmap``,
    %           otherwise, it is ``false``.
    %           (**optional**. If missing, ``val`` will remain an empty array ``[]`` on output.)
    %
    %   Interface
    %   ---------
    %
    %       [val, failed] = pm.matlab.hashmap.getVal(key, hashmap)
    %
    %   Example
    %   -------
    %
    %       hashmap = {"key1", 1, "key2", "val2", "key3", false};
    %       val = pm.matlab.hashmap.getVal("key2", hashmap) % = "val2"
    %       val = pm.matlab.hashmap.getVal("key3", hashmap) % = false
    %       val = pm.matlab.hashmap.getVal("key3", hashmap(1:4)) % = {}
    %       val = pm.matlab.hashmap.getVal("key2", hashmap(1:4)) % = "val2"
    %       val = pm.matlab.hashmap.getVal("key2", hashmap(1:3)) % error
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    vararginLen = length(hashmap);
    if mod(vararginLen, 2) ~= 0
        help("pm.matlab.hashmap.getVal");
        error   ( newline ...
                + "The length of ``hashmap`` must be even." + newline ...
                + newline ...
                + "length(hashmap) = " + string(vararginLen) + newline ...
                + newline ...
                );
    end
    val = [];
    if 1 < nargout
        failed = true;
    end
    for i = 1 : vararginLen - 1
        if strcmpi(string(hashmap{i}), string(key))
            if 1 < nargout
                failed = false;
            end
            val = hashmap{i + 1};
            return;
        end
    end
end
function [val, keyFound] = getKeyVal(key,varargin)
    vararginLen = length(varargin);
    keyString = lower(string(key));
    val = {};
    keyFound = false;
    for i = 1:2:vararginLen-1
        if strcmp(lower(string(varargin{i})),keyString)
            keyFound = true;
            val = varargin{i+1};
            break;
        end
    end
end
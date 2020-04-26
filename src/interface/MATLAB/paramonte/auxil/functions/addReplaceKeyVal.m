function newKeyWordSet = addReplaceKeyVal(key,val,varargin)
    newKeyWordSet = {varargin{:}};
    vararginLen = length(newKeyWordSet);
    keyString = string(key);
    keyFound = false;
    for i = 1:vararginLen-1
        try
            if strcmpi(string(newKeyWordSet{i}),keyString)
                keyFound = true;
                newKeyWordSet{i+1} = val;
                break;
            end
        catch
            continue
        end
    end
    if ~keyFound
        newKeyWordSet = { newKeyWordSet{:} , key, val };
    end
end
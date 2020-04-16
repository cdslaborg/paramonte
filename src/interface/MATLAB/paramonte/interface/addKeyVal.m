function newKeyWordSet = addKeyVal(key,val,varargin)
    newKeyWordSet = {varargin{:}};
    [currentVal, keyFound] = getKeyVal(key,varargin{:});
    if ~keyFound
        newKeyWordSet = [ newKeyWordSet , {key,val} ];
    end
end
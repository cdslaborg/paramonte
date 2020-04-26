function newKeyValSet = addKeyVal(key,val,varargin)
    newKeyValSet = varargin;
    if nargin>2 && getVecLen(newKeyValSet{1})
        [currentVal, keyFound] = getKeyVal(key,varargin{:});
        if ~keyFound
            newKeyValSet = { newKeyValSet{:} , key, val };
        end
    else
        newKeyValSet = { key, val };
    end
end
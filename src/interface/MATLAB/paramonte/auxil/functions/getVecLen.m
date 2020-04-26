function len = getVecLen(object)
    % char vector has a length of zero or one
    % empty string has a length of zero
    len = 0;
    if isa(object,"char"); object = string(object); end
    for i = 1:length(object)
        try 
            if ~strcmp(string(object(i)),"")
                len = len + 1;
            end
        catch
            len = len + 1;
        end
    end
end
function enclosedString = encloseString(inputString)
    % enclose string with either single or double quotation. 
    % Return None if string contains both quotation symbols
    hasSingleQuote = contains(inputString,"'");
    hasDoubleQuote = contains(inputString,'"');
    if hasSingleQuote
        if hasDoubleQuote
            enclosedString = [];
        else
            enclosedString = '"' + inputString + '"';
        end
    else
        enclosedString = "'" + inputString + "'";
    end
end

function result = isNumericString(object)
    result = ~isnan(str2double(object));
end
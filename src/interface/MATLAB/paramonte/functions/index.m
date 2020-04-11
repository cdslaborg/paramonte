function index_val = index(str, sub_str)

    index_val = strfind(str, sub_str);
    
    if isempty(index_val)
        index_val = 0; 
    else
        index_val = index_val(1); 
    end
    
end
function cellArray = convertStruct2Cell(structure,excludeList)
    cellArray = {};
    fnameList = fieldnames(structure);
    excludeList = string(excludeList);
    for i = 1:length(fnameList)
        fname = string(fnameList{i});
        if ~any(strcmp(excludeList,fname)) 
            if ~isempty(structure.(fname))
                cellArray = { cellArray{:}, fname, structure.(fname) };
            end
        end
    end
end
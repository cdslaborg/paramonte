function newto = copyComponents(from,to,fromFieldList)
    %
    % Copy the contents of the struct/object ``from`` to the struct/object ``to`` 
    % recursively and without destroying the existing components in ``to``.
    %
    if isempty(to); to = struct(); end
    if nargin<3
        try
            fromFieldList = fieldnames(from);
        catch
            try
                fromFieldList = properties(from);
            catch
                fromFieldList = [];
            end
        end
    end

    if isempty(fromFieldList)
        newto = from;
        return;
    else
        newto = to;
    end

    for i = 1:length(fromFieldList)

        % check if the field is itself a struct or object

        try
            subFieldList = fieldnames(from.(fromFieldList{i}));
        catch
            try
                subFieldList = properties(from.(fromFieldList{i}));
            catch
                subFieldList = [];
            end
        end

        if isempty(subFieldList)
            newtoComponent = from.(fromFieldList{i});
        else
            try
                newtoComponent = newto.(fromFieldList{i});
            catch % the field of from does not exist in to.
                newtoComponent = struct();
            end
            newtoComponent = copyComponents(from.(fromFieldList{i}), newtoComponent, subFieldList);
        end

        try % struct
            newto.(fromFieldList{i}) = newtoComponent;
        catch % object
            newto.addprop(fromFieldList{i});
            newto.(fromFieldList{i}) = newtoComponent;
        end

end
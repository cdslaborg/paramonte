function tonew = copy(from, to, fields)
    %
    %   Copy the contents of the struct/object ``from``
    %   to the struct/object ``to`` recursively and without
    %   destroying the existing components in ``to``.
    %
    %   Parameters
    %   ----------
    %
    %       from
    %
    %           The input scalar MATLAB struct whose (select)
    %           components must be copy/merged with the components
    %           of the input struct ``to``.
    %
    %       to
    %
    %           The input scalar MATLAB struct to which the
    %           components of ``from`` struct must be copied.
    %           (**optional**, default = ``struct()``)
    %
    %       fields
    %
    %           The input vector of MATLAB strings each element
    %           of which is a field name in ``from`` whose value
    %           has to be copied to the struct ``to``.
    %           (**optional**, default = ``fieldnames(from)`` or ``properties(from)``)
    %
    %   Returns
    %   -------
    %
    %       tonew
    %
    %           The output MATLAB struct containing the
    %           merger of the two input MATLAB structs.
    %           If a field name in ``fields`` is common between
    %           ``from`` and ``to``, the field value of ``from`` will
    %           overwrite the corresponding field value of ``to``
    %           in the output ``tonew``.
    %
    %   Interface
    %   ---------
    %
    %       tonew = pm.matlab.struct.copy(from)
    %       tonew = pm.matlab.struct.copy(from, to)
    %       tonew = pm.matlab.struct.copy(from, [], fields)
    %       tonew = pm.matlab.struct.copy(from, to, fields)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if nargin < 3
        try
            fields = fieldnames(from);
        catch
            try
                fields = properties(from);
            catch
                fields = [];
            end
        end
    end
    if nargin < 2
        to = [];
    end
    if isempty(to)
        to = struct();
    end
    if isempty(fields)
        tonew = from;
        return;
    else
        tonew = to;
    end

    for i = 1 : length(fields)
        % check if the ``field`` is itself a struct or object.
        try
            subFieldList = fieldnames(from.(fields{i}));
        catch
            try
                subFieldList = properties(from.(fields{i}));
            catch
                subFieldList = [];
            end
        end
        if isempty(subFieldList)
            tonewComponent = from.(fields{i});
        else
            try
                tonewComponent = tonew.(fields{i});
            catch % the ``field`` of from does not exist in to.
                tonewComponent = struct();
            end
            tonewComponent = pm.matlab.struct.copy(from.(fields{i}), tonewComponent, subFieldList);
        end
        try % struct
            tonew.(fields{i}) = tonewComponent;
        catch % object
            tonew.addprop(fields{i});
            tonew.(fields{i}) = tonewComponent;
        end
    end
end
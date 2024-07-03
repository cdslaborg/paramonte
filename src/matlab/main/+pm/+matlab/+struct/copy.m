%>  \brief
%>  Copy the contents of the struct/object ``from``
%>  to the struct/object ``to`` recursively and without
%>  destroying the existing components in ``to``.<br>
%>
%>  \param[in]  from    :   The input scalar MATLAB struct whose (select)
%>                          components must be copy/merged with the components
%>                          of the input struct ``to``.<br>
%>  \param[in]  to      :   The input scalar MATLAB struct to which the
%>                          components of ``from`` struct must be copied.<br>
%>                          (**optional**, default = ``struct()``)
%>  \param[in]  fields  :   The input vector of MATLAB strings each element
%>                          of which is a field name in ``from`` whose value
%>                          has to be copied to the struct ``to``.<br>
%>                          (**optional**, default = ``fieldnames(from)`` or ``properties(from)``)
%>
%>  \return
%>  ``tonew``           :   The output MATLAB struct containing the
%>                          merger of the two input MATLAB structs.<br>
%>                          If a field name in ``fields`` is common between
%>                          ``from`` and ``to``, the field value of ``from`` will
%>                          overwrite the corresponding field value of ``to``
%>                          in the output ``tonew``.<br>
%>
%>  \interface{copy}
%>  \code{.m}
%>
%>      tonew = pm.matlab.struct.copy(from)
%>      tonew = pm.matlab.struct.copy(from, to)
%>      tonew = pm.matlab.struct.copy(from, [], fields)
%>      tonew = pm.matlab.struct.copy(from, to, fields)
%>
%>  \endcode
%>
%>  \example{copy}
%>  \include{lineno} example/matlab/struct/copy/main.m
%>  \output{copy}
%>  \include{lineno} example/matlab/struct/copy/main.out.m
%>
%>  \final{copy}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:59 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function tonew = copy(from, to, fields)
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
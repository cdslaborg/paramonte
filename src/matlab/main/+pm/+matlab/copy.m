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
%>  \param[in]  field   :   The input vector of MATLAB strings each element
%>                          of which is a field or property name in ``from``
%>                          whose value has to be copied to the struct ``to``.<br>
%>                          (**optional**, default = ``fieldnames(from)`` or ``properties(from)``)
%>  \param[in]  exclude :   The input vector of MATLAB strings each element
%>                          of which is a field or property name in ``from``
%>                          whose value has to be skipped and excluded from the copy to the struct ``to``.<br>
%>                          (**optional**, default = ``[]``)
%>
%>  \return
%>  ``tonew``           :   The output MATLAB struct containing the
%>                          merger of the two input MATLAB structs.<br>
%>                          If a field name in ``field`` is common between
%>                          ``from`` and ``to``, the field value of ``from`` will
%>                          overwrite the corresponding field value of ``to``
%>                          in the output ``tonew``.<br>
%>
%>  \interface{copy}
%>  \code{.m}
%>
%>      tonew = pm.matlab.copy(from)
%>      tonew = pm.matlab.copy(from, to)
%>      tonew = pm.matlab.copy(from, [], field)
%>      tonew = pm.matlab.copy(from, [], field, [])
%>      tonew = pm.matlab.copy(from, to, field, [])
%>      tonew = pm.matlab.copy(from, to, field, exclude)
%>      tonew = pm.matlab.copy(from, to, [], exclude)
%>      tonew = pm.matlab.copy(from, [], [], exclude)
%>
%>  \endcode
%>
%>  \warning
%>  If a name appears in both ``field`` and ``exclude`` input arguments, then it is excluded.<br>
%>
%>  \example{copy}
%>  \include{lineno} example/matlab/copy/main.m
%>  \output{copy}
%>  \include{lineno} example/matlab/copy/main.out.m
%>
%>  \final{copy}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:59 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function tonew = copy(from, to, field, exclude)

    if  nargin < 4
        exclude = [];
    elseif pm.array.len(exclude) == 0
        exclude = [];
    else
        exclude = [string(exclude)];
    end

    if  nargin < 3
        field = [];
    end
    if  isempty(field)
        try
            field = fieldnames(from);
        catch
            try
                field = properties(from);
            catch
                field = [];
            end
        end
    end

    if  nargin < 2
        to = [];
    end
    if  isempty(to)
        to = struct();
    end

    if isempty(field)
        tonew = from;
        return;
    else
        tonew = to;
    end

    for i = 1 : length(field)

        %%%%
        %%%% Skip the ``field`` if it is excluded.
        %%%%

        if  isempty(exclude) || ~any(strcmp(exclude, field(i)))

            %%%%
            %%%% Check if the ``field`` is itself a struct or object.
            %%%%

            try
                subFieldList = fieldnames(from.(field{i}));
            catch
                try
                    subFieldList = properties(from.(field{i}));
                catch
                    subFieldList = [];
                end
            end

            %%%%
            %%%% Create the new copy component.
            %%%%

            if  isempty(subFieldList)
                tonewComponent = from.(field{i});
            else
                try
                    tonewComponent = tonew.(field{i});
                catch % the ``field`` of from does not exist in to.
                    tonewComponent = struct();
                end
                tonewComponent = pm.matlab.copy(from.(field{i}), tonewComponent, subFieldList);
            end

            %%%%
            %%%% Copy the field values.
            %%%%

            try % struct
                tonew.(field{i}) = tonewComponent;
            catch % object
                tonew.addprop(field{i});
                tonew.(field{i}) = tonewComponent;
            end

        end

    end

end
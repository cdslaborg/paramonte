%>  \brief
%>  Return a copy of the input ``object``
%>  whose property (or field) values are overwritten
%>  with the corresponding ``(key, val)`` pairs in the
%>  input ``hashmap`` cell array.
%>
%>  \details
%>  This functionality is useful for converting
%>  the variadic arguments of functions to the
%>  corresponding property values of objects
%>  with the functions.
%>
%>  \param[in]  hashmap     :   The input cell array of even number of elements
%>                              containing the field names and values of the input ``object``
%>                              as ``(key, val)`` pairs stored sequentially as the cell elements.<br>
%>                              All keys in the input ``hashmap`` must already exist as
%>                              fields/properties of the input object.<br>
%>  \param[in]  object      :   The input scalar MATLAB struct
%>                              or object whose property/field values
%>                              are to be overwritten with the corresponding
%>                              values from the input ``hashmap`` cell.<br>
%>  \param[in]  insensitive :   The input scalar MATLAB logical.
%>                              If ``false``, then keys within the input ``hashmap`` will be matched
%>                              against the input ``object`` properties(or fields) case-sensitively.<br>
%>                              (**optional**, default = ``false``)
%>  \param[in]  extensible  :   The input scalar MATLAB logical.<br>
%>                              If ``true``, then keys in the input ``hashmap`` that are missing in the
%>                              input ``object`` properties(or fields) will be added to output ``objnew``.<br>
%>                              This functionality requires the input ``object`` to be either a MATLAB ``struct``
%>                              or an object whose ultimate superclass is the MATLAB ``handle`` class,
%>                              which allows adding new properties via ``addprop()`` method.<br>
%>                              (**optional**, default = ``false``)
%>  \param[in]  recursive   :   The input scalar MATLAB logical.<br>
%>                              If ``true``, then key-val pairs in the input ``hashmap`` that match struct-cell
%>                              or match object-cell pattern will be also converted recursively like the parent
%>                              input ``object`` until all sub-components are obtained recursively.<br>
%>                              (**optional**, default = ``false``)
%>
%>  \return
%>  objnew                  :   The output copy of the input object whose
%>                              property/field values are overwritten with the
%>                              corresponding values from the input ``hashmap`` cell.<br>
%>
%>  \interface{hash2comp}
%>  \code{.m}
%>
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object)
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, [])
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, [], [])
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, [], [], [])
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, [], [], recursive)
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, [], extensible, [])
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, insensitive, [], [])
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, [], extensible, recursive)
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, insensitive, [], recursive)
%>      objnew = pm.matlab.hashmap.hash2comp(hashmap, object, insensitive, extensible, recursive)
%>
%>  \endcode
%>
%>  \example{hash2comp}
%>  \include{lineno} example/matlab/hashmap/hash2comp/main.m
%>  \output{hash2comp}
%>  \include{lineno} example/matlab/hashmap/hash2comp/main.out.m
%>
%>  \todo
%>  \pvhigh
%>  Is the last condition ``~isempty(properties(objnew.(key)))`` in the implementation of this procedure necessary?<br>
%>  Its presence contradicts the function behavior at the zeroth level of iteration and subsequent iterations.<br>
%>  This issue must be addressed as soon as the correction is verified to not change the behavior of other
%>  functions calling this function.<br>
%>
%>  \final{hash2comp}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:47 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function objnew = hash2comp(hashmap, object, insensitive, extensible, recursive)
    objnew = object;
    hashmapLen = length(hashmap);
    if  rem(hashmapLen, 2) ~= 0
        help("pm.matlab.hashmap.hash2comp");
        disp("length(hashmap)");
        disp( length(hashmap) );
        error   ( newline ...
                + "The input argument ``hashmap`` must be a cell array of even length." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end
    if ~isstruct(objnew)
        fieldList = properties(objnew);
        fieldListLen = length(fieldList);
    else
        fieldList = fieldnames(objnew);
        fieldListLen = length(fieldList);
    end
    if  nargin < 5
        recursive = [];
    end
    if  nargin < 4
        extensible = [];
    end
    if  nargin < 3
        insensitive = [];
    end
    if  isempty(recursive)
        recursive = false;
    end
    if  isempty(extensible)
        extensible = false;
    end
    if  isempty(insensitive)
        insensitive = false;
    end
    for i = 1 : 2 : hashmapLen
        keyfound = false;
        key = string(hashmap{i});
        for ip = 1 : fieldListLen
            if ~insensitive
                keyfound = strcmp(key, string(fieldList(ip)));
            else
                keyfound = strcmpi(key, string(fieldList(ip)));
            end
            if  keyfound
                key = fieldList{ip};
                break;
            end
        end
        if ~keyfound
            if ~extensible
                error   ( newline ...
                        + "The requested ``object`` property:" + newline ...
                        + newline ...
                        +  pm.io.tab() + """" + string(hashmap{i}) + """" + newline ...
                        + newline ...
                        + "does not exist in the specified input ``object``." + newline ...
                        + newline ...
                        );
            elseif ~isstruct(objnew)
                objnew.addprop(key);
            end
        end
        % \todo:
        % Is the last condition ``~isempty(properties(objnew.(key)))`` necessary?
        % Its presence contradicts the function behavior at the zeroth level of iteration and subsequent iterations.
        % This issue must be addressed as soon as it is clear the correction does not change the behavior of other
        % functions calling this function.
        if  recursive && isa(hashmap{i + 1}, "cell") && (isa(objnew.(key), "struct") || isa(objnew.(key), "handle") || ~isempty(properties(objnew.(key))))
            objnew.(key) = pm.matlab.hashmap.hash2comp(hashmap{i + 1}, objnew.(key), insensitive, extensible, recursive);
        else
            objnew.(key) = hashmap{i + 1};
        end
    end
end
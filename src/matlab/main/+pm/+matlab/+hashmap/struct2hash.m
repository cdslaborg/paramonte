%>  \brief
%>  Return a ``hashmap`` cell array containing all field names
%>  and field values of input scalar MATLAB ``object`` as ``(key, val)``
%>  stored sequentially in the cell array.<br>
%>
%>  \param[in]  object      :   The input scalar MATLAB struct
%>                              whose field names and values are
%>                              to be converted to a ``hashmap`` cell.<br>
%>  \param[in]  exkeys      :   The input vector of MATLAB strings
%>                              containing a list of field names of the input ``object``
%>                              that must be excluded from the output ``hashmap`` cell.<br>
%>                              If the input argument ``unique`` is set to ``True``,
%>                              then all elements of ``exkeys`` are compared to the
%>                              object field and subfield names case-insensitively.<br>
%>                              (**optional**, default = ``[]``)
%>  \param[in]  unique      :   The input scalar MATLAB logical.<br>
%>                              If ``true``, only the first instance of occurrence
%>                              of any field name is added to the output cell array
%>                              without considering case-sensitivity of the field names.<br>
%>                              This argument also affects the input argument ``exkeys``.<br>
%>                              (**optional**, default = ``false``)
%>  \param[in]  onlyfull    :   The input scalar MATLAB logical.<br>
%>                              If ``true``, only the field names with field values that
%>                              are nonempty (as assessed by the MATLAB intrinsic ``isempty()``)
%>                              are included in the output ``hashmap``.<br>
%>                              (**optional**, default = ``false``)
%>
%>  \return
%>  ``hashmap``             :   The output cell array of even number of elements
%>                              containing the field names and values of the input ``object``
%>                              as ``(key, val)`` pairs stored sequentially as the cell elements.<br>
%>
%>  \interface{struct2hash}
%>  \code{.m}
%>
%>      hashmap = pm.matlab.hashmap.struct2hash(object)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, exkeys)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, [], unique)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, exkeys, unique)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, exkeys, [], onlyfull)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, [], unique, onlyfull)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, exkeys, unique, onlyfull)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, [], [], onlyfull)
%>      hashmap = pm.matlab.hashmap.struct2hash(object, [], [], [])
%>
%>  \endcode
%>
%>  \example{struct2hash}
%>  \include{lineno} example/matlab/hashmap/struct2hash/main.m
%>  \output{struct2hash}
%>  \include{lineno} example/matlab/hashmap/struct2hash/main.out.m
%>
%>  \final{struct2hash}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:56 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hashmap = struct2hash(object, exkeys, unique, onlyfull)
    fnameList = fieldnames(object);
    fnameListLen = length(fnameList);
    if nargin < 4; onlyfull = []; end
    if nargin < 3; unique = []; end
    if nargin < 2; exkeys = []; end
    if isempty(unique); unique = false; end
    if isempty(onlyfull); onlyfull = false; end
    hashmap = cell(fnameListLen * 2, 1);
    exkeys = string(exkeys);
    if  unique
        fnameListLower = lower(fnameList);
        exkeys = lower(exkeys);
    end
    counter = 0;
    for i = 1 : fnameListLen
        fname = string(fnameList{i});
        if ~any(strcmp(exkeys, fname))
            if ~(onlyfull && isempty(object.(fname))) && ~(unique && any(strcmpi(fnameListLower(1 : i - 1), fname)))
                counter = counter + 1;
                hashmap{counter} = fname;
                counter = counter + 1;
                hashmap{counter} = object.(fname);
            end
        end
    end
    hashmap = hashmap(1 : counter);
end
%>  \brief
%>  Return the input ``hashmap`` with the new ``(key, val)`` appended,
%>  only if the input pair does not exist in the input ``hashmap``.<br>
%>  If the pair exists, replace it with the new value.<br>
%>
%>  \param[in]  key     :   The input scalar MATLAB string
%>                          containing the key to add to the input pair list.
%>  \param[in]  val     :   The input value to be appended to the
%>                          input ``hashmap`` after the input ``key``.
%>  \param[in]  hashmap :   The input cell array of even number of elements
%>                          containing the ``(key, val)`` pairs of the ``hashmap``
%>                          in sequence as element of the cell array.
%>
%>  \return
%>  ``hashnew``         :   The output cell array of even number of elements
%>                          containing the input pair list and ``(key, val)`` pair.<br>
%>                          If the input ``key`` exists in the input ``hashmap``,
%>                          the input ``(key, val)`` pair will not be added.<br>
%>
%>  \interface{repKeyVal}
%>  \code{.m}
%>
%>      hashnew = pm.matlab.hashmap.repKeyVal(key, val, hashmap)
%>
%>  \endcode
%>
%>  \example{repKeyVal}
%>  \include{lineno} example/matlab/hashmap/repKeyVal/main.m
%>  \output{repKeyVal}
%>  \include{lineno} example/matlab/hashmap/repKeyVal/main.out.m
%>
%>  \final{repKeyVal}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:52 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hashnew = repKeyVal(key, val, hashmap)
    hashnew = hashmap;
    vararginLen = length(hashnew);
    keyString = string(key);
    failed = true;
    for i = 1 : vararginLen - 1
        try
            if strcmpi(string(hashnew{i}), keyString)
                hashnew{i + 1} = val;
                failed = false;
                break;
            end
        catch
            continue
        end
    end
    if  failed
        hashnew = {hashnew{:}, key, val};
    end
end
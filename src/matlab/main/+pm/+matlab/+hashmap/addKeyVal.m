%>  \brief
%>  Return the input hashmap with the new (key, val) appended,
%>  only if the input pair does not exist in the input hashmap.
%>
%>  \param[in]  key     :   The input scalar MATLAB string
%>                          containing the key to add to the input pair list.
%>  
%>  \param[in]  val     :   The input value to be appended to the 
%>                          input hashmap after the input ``key``.
%>  
%>  \param[in]  hashmap :   The input cell array of even number of elements
%>                          containing a set of ``(key, val)`` pairs of the hashmap.
%>  
%>  \return
%>  `hashnew`           :   The output cell array of even number of elements
%>                          containing the input pair list and ``(key, val)`` pair.
%>                          If the input ``key`` exists in the input ``hashmap``,
%>                          the input ``(key, val)`` pair will not be added.
%>
%>  \interface{addKeyVal}
%>  \code{.m}
%>
%>      hashnew = pm.matlab.hashmap.addKeyVal(key, val, hashmap)
%>
%>  \endcode
%>  \example{addKeyVal}
%>
%>      hashmap = {"key1", 1, "key2", "val2"};
%>      hashnew = pm.matlab.hashmap.addKeyVal("key3", false, {})
%>      hashnew = pm.matlab.hashmap.addKeyVal("key2", "val2", hashmap) % = hashmap
%>      hashnew = pm.matlab.hashmap.addKeyVal("key3", false, hashmap)
%>
%>  \final{addKeyVal}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:39 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function hashnew = addKeyVal(key, val, hashmap)
    hashnew = hashmap;
    if  2 < nargin && 0 < length(hashnew)
        [currentVal, failed] = pm.matlab.hashmap.getKeyVal(key, hashmap);
        if  failed
            hashnew = {hashnew{:}, key, val};
        end
    else
        hashnew = {key, val}; % XXX is this correct behavior?
    end
end
%>  \brief
%>  Return the value corresponding to the input ``key``
%>  in the input ``hashmap`` cell array.
%>
%>  \param[in]  key     :   The input scalar MATLAB string
%>                          containing the key to search for in the input pair list.
%>
%>  \param[in]  hashmap :   The input cell array of even number of elements
%>                          containing the ``(key, val)`` pairs of the hashmap
%>                          in sequence as element of the cell array.
%>
%>  \return
%>  ``val``             :   The output MATLAB object stored in the element of the
%>                          input cell array, whose ``key`` is given as the input.
%>                          <br>
%>  ``failed``          :   The output scalar MATLAB logical that is ``true``
%>                          if and only if the input ``key`` exists in the input ``hashmap``,
%>                          otherwise, it is ``false``.<br>
%>                          (**optional**. If missing, ``val`` will remain an empty array ``[]`` on output.)
%>
%>  \interface{getKeyVal}
%>  \code{.m}
%>
%>      val = pm.matlab.hashmap.getKeyVal(key, hashmap)
%>      [val, failed] = pm.matlab.hashmap.getKeyVal(key, hashmap)
%>
%>  \endcode
%>
%>  \warning
%>  The length of the input ``hashmap`` must be even.<br>
%>  Ann odd length of the input cell array will lead to a runtime call to the MATLAB intrinsic ``error()``.<br>
%>
%>  \example{getKeyVal}
%>  \include{lineno} example/matlab/hashmap/getKeyVal/main.m
%>  \output{getKeyVal}
%>  \include{lineno} example/matlab/hashmap/getKeyVal/main.out.m
%>
%>  \final{getKeyVal}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:43 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [val, failed] = getKeyVal(key, hashmap)
    vararginLen = length(hashmap);
    if mod(vararginLen, 2) ~= 0
        help("pm.matlab.hashmap.getKeyVal");
        error   ( newline ...
                + "The length of ``hashmap`` must be even." + newline ...
                + newline ...
                + "length(hashmap) = " + string(vararginLen) + newline ...
                + newline ...
                );
    end
    val = [];
    if  1 < nargout
        failed = true;
    end
    for i = 1 : vararginLen - 1
        if strcmpi(string(hashmap{i}), string(key))
            if 1 < nargout
                failed = false;
            end
            val = hashmap{i + 1};
            return;
        end
    end
end
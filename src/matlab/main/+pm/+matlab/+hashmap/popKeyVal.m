%>  \brief
%>  Return the input ``hashmap`` while all cases
%>  of keys matching the input ``keys`` are deleted
%>  along with their corresponding values immediately
%>  appearing after the elements matching the input ``keys``.
%>
%>  \param[in]  keys    :   The input MATLAB cell array of rank 1 containing
%>                          the key(s) to search for in the input pair list.
%>  \param[in]  hashmap :   The input cell array of even number of elements
%>                          containing the ``(key, val)`` pairs of the input
%>                          ``hashmap`` in sequence as element of the cell array.
%>  
%>  \return
%>  ``keyval``          :   The output cell array of even number of elements
%>                          containing the ``(key, val)`` pairs of the keys
%>                          in the input argument ``keys``.<br>
%>  ``hashout``         :   The output cell array of even number of elements
%>                          containing the ``(key, val)`` pairs of the input
%>                          ``hashmap`` while all instances of keys matching
%>                          the input ``keys`` along with their values are deleted.<br>
%>
%>  \interface{popKeyVal}
%>  \code{.m}
%>
%>      [keyval, hashout] = pm.matlab.hashmap.popKeyVal(keys, hashmap)
%>
%>  \endcode
%>
%>  \example{popKeyVal}
%>  \include{lineno} example/matlab/hashmap/popKeyVal/main.m
%>  \output{popKeyVal}
%>  \include{lineno} example/matlab/hashmap/popKeyVal/main.out.m
%>
%>  \final{popKeyVal}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 10:49 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [keyval, hashout] = popKeyVal(keys, hashmap)
    vararginLen = length(hashmap);
    if mod(vararginLen, 2) ~= 0
        help("pm.matlab.hashmap.popKeyVal");
        error   ( newline ...
                + "The length of ``hashmap`` must be even." + newline ...
                + newline ...
                + "length(hashmap) = " + string(vararginLen) + newline ...
                + newline ...
                );
    end
    hashout = cell(length(hashmap), 1);
    keyval = cell(length(hashmap), 1);
    counter = 0;
    kvcount = 0;
    strkeys = string(keys);
    for jkey = 1 : 2 : vararginLen
        if ~any(strcmpi(string(hashmap{jkey}), strkeys))
            counter = counter + 1;
            hashout{counter} = hashmap{jkey};
            counter = counter + 1;
            hashout{counter} = hashmap{jkey + 1};
        else
            kvcount = kvcount + 1;
            keyval{kvcount} = hashmap{jkey};
            kvcount = kvcount + 1;
            keyval{kvcount} = hashmap{jkey + 1};
        end
    end
    hashout = hashout(1 : counter);
    keyval = keyval(1 : kvcount);
end
%>  \brief
%>  Return the lists of values and indices of elements
%>  of the input ``strlist`` that match the input scalar
%>  (or vector of) string(s) or integer(s) index ``keylist``.<br>
%>
%>  \details
%>  Beware that all specified values in the input ``keylist``
%>  must exist in the input argument ``strlist``, otherwise,
%>  the procedure will interrupt the program by raising an exception.
%>
%>  \param[in]  strlist :   The input scalar (or vector of) MATLAB string(s).<br>
%>  \param[in]  keylist :   The input scalar (or vector of) MATLAB string(s)
%>                          or cell array of char vectors or vector of MATLAB
%>                          integers or a cell array of mix of the element
%>                          types to match the input ``strlist``.<br>
%>
%>  \return
%>  ``loclist``         :   The output scalar (or vector of same size as ``keylist`` of)
%>                          MATLAB integer(s) containing the location(s) of the occurrence(s)
%>                          of the input ``keylist``.<br>
%>                          <br>
%>  ``namlist``         :   The output scalar (or vector of same size as ``keylist`` of)
%>                          MATLAB string(s) containing the name(s) of the occurrence(s)
%>                          of the input ``keylist``.<br>
%>
%>  \interface{locname}
%>  \code{.m}
%>
%>      [loclist, namlist] = pm.str.locname(strlist, keylist)
%>
%>  \endcode
%>
%>  \example{locname}
%>  \include{lineno} example/str/locname/main.m
%>  \output{locname}
%>  \include{lineno} example/str/locname/main.out.m
%>
%>  \final{locname}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:38 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [loclist, namlist] = locname(strlist, keylist)

    keylistLen = pm.array.len(keylist);
    strlist = string(strlist);
    if keylistLen == 0
        loclist = 1 : 1 : length(strlist);
        namlist = strlist;
        return;
    end

    failed = false;
    strlistLen = length(strlist);
    for i = keylistLen : -1 : 1
        failed = true;
        key = keylist(i);
        if  isa(key, "cell")
            key = key{1};
        end
        if isa(key, "string") || isa(key, "char")
            %namlist(i) = key;
            for j = 1 : strlistLen
                if strcmpi(key, strlist(j))
                    namlist(i) = strlist(j);
                    loclist(i) = j;
                    failed = false;
                    break;
                end
            end
            if failed
                % This situation happens because of automatic conversion of numbers
                % to strings in MATLAB array constructs like a = ["paramonte", 1];
                key = str2double(key);
                failed = isnan(key);
                if ~failed
                    failed = key < 1 || strlistLen < key;
                    if ~failed
                        loclist(i) = key;
                        namlist(i) = strlist(key);
                    end
                end
            end
        elseif isnumeric(key)
            loclist(i) = key;
            failed = strlistLen < loclist(i);
            if ~failed
                namlist(i) = string(strlist{loclist(i)});
            end
        end
        if  failed
            break;
        end
    end

    if  failed
        help("pm.str.locname");
        disp("strlist = ");
        disp(strlist);
        disp("keylist = ");
        disp(keylist);
        error   ( newline ...
                + "The element #" + string(i) + " of the input ``keylist`` above" + newline ...
                + "argument does not match any elements of the input ``strlist`` above." + newline ...
                + newline ...
                );
    end

end

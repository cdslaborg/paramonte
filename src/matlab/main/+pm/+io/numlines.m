%>  \brief
%>  Return a scalar MATLAB integer representing
%>  the number of lines in the specified input ``idname``.
%>
%>  \param[in]  idname  :   The input scalar MATLAB object that can be either:<br>
%>                          <ol>
%>                              <li>    a MATLAB string containing the path
%>                                      to a record-based external file whose
%>                                      number of records (lines) must be returned.<br>
%>                              <li>    a MATLAB integer, representing the file
%>                                      ID of an opened record-based file whose
%>                                      number of records (lines) must be returned.<br>
%>                          </ol>
%>
%>  \return
%>  ``nlines``          :   The output scalar MATLAB integer containing the
%>                          number of records (lines) in the specified input ``idname``.
%>
%>  \interface{numlines}
%>  \code{.m}
%>
%>      nlines = pm.io.numlines(idname);
%>
%>  \endcode
%>
%>  \warning
%>  The input file (ID) will be closed even if it is open upon entry to the function.
%>
%>  \example{numlines}
%>  \include{lineno} example/io/numlines/main.m
%>  \output{numlines}
%>  \include{lineno} example/io/numlines/main.out.m
%>
%>  \final{numlines}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:50 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function nlines = numlines(idname)
    if  pm.introspection.verified(idname, "string", 1)
        [fid, errmsg] = fopen(idname, 'r');
        if fid < 0
            help("pm.io.numlines");
            error( newline + ...
                + "Failed to open file """ + idname + """: " + newline ...
                + newline ...
                + string(errmsg) + newline ...
                + newline ...
                + "For more information, see the documentation of the function displayed above." + newline ...
                + newline ...
                );
        end
    elseif pm.introspection.verified(idname, "integer", 1)
        fid = idname;
    else
        help("pm.io.numlines");
        disp("idname = ");
        disp(idname);
        error   ( newline ...
                + "Invalid input value for ``idname`` displayed above." + newline ...
                + "For more information, see the documentation of the function displayed above." + newline ...
                + newline ...
                );
    end
    nlines = 0;
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        else
            nlines = nlines + 1;
        end
    end
    fclose(fid);
end
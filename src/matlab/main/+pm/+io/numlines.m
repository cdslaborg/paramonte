function nlines = numlines(idname)
    %
    %   Return a scalar MATLAB integer representing
    %   the number of lines in the specified input ``idname``.
    %
    %   Parameters
    %   ----------
    %
    %       idname
    %
    %           The input scalar MATLAB object that can be either:
    %
    %               1.  a MATLAB string containing the path
    %                   to a record-based external file whose
    %                   number of records (lines) must be returned.
    %
    %               2.  a MATLAB integer, representing the file
    %                   ID of an opened record-based file whose
    %                   number of records (lines) must be returned.
    %
    %   Returns
    %   -------
    %
    %       nlines
    %
    %           The output scalar MATLAB integer containing the
    %           number of records (lines) in the specified input ``idname``.
    %
    %   Interface
    %   ---------
    %
    %       nlines = pm.io.numlines(idname);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if pm.introspection.istype(idname, "string", 1)
        [fid, errmsg] = fopen(idname, 'r');
        if fid < 0
            help("pm.io.numlines");
            error( newline + ...
                + "Failed to open file """ + idname + """: " + newline ...
                + newline ...
                + string(errmsg) + newline ...
                + newline ...
                + "For more information, see the documentation of the function displayed above" + newline ...
                + newline ...
                );
        end
    elseif pm.introspection.istype(idname, "integer", 1)
        fid = idname;
    else
        help("pm.io.numlines");
        disp("idname = ");
        disp(idname);
        error   ( newline ...
                + "Invalid input value for ``idname`` displayed above." + newline ...
                + "For more information, see the documentation of the function displayed above" + newline ...
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
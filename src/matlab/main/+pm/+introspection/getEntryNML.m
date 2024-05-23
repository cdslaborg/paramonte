%>  \brief
%>  Return a Fortran-namelist-compatible conversion of the input ``varval``.
%>
%>  \details
%>  This functionality is primarily used by the ParaMonte MATLAB internal
%>  routines to communicate information with Fortran shared libraries.
%>  As such, it is of limited to most end users of the library.
%>
%>  \param[in]  varname :   The input scalar MATLAB string containing the label to assign
%>                          to the namelist-converted value in the output ``entry``.<br>
%>                          The specified value of ``varname`` will be trimmed
%>                          (to remove leading and trailing blanks).
%>
%>  \param[in]  varval  :   The input value to be converted to namelist-compatible value.
%>
%>  \param[in]  vartype :   See the documentation of the corresponding
%>                          argument of ``pm.introspection.istype()``.
%>
%>  \param[in]  varsize :   See the documentation of the corresponding
%>                          argument of ``pm.introspection.istype()``.
%>
%>  \return
%>  `entry`             :   The output scalar MATLAB string containing the namelist-compatible
%>                          conversion of the input value ``varval`` and the given ``varname``
%>                          in the following format:    varname=namelist-compatible-varval
%>
%>                          \note
%>                          If the input ``varval`` is an array, its elements will be comma-separated.
%>                          If the input value is string, it will be quoted properly.
%>
%>  \interface{getEntryNML}
%>  \code{.m}
%>
%>      entry = pm.introspection.getEntryNML(varname, varval, vartype, varsize)
%>
%>  \endcode
%>  \final{getEntryNML}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:38 PM, University of Texas at Arlington<br>
%>
function entry = getEntryNML(varname, varval, vartype, varsize)
    varval = varval(:);
    varvalen = numel(varval);
    %varname = string(inputname(2));
    varname = string(strtrim(varname));
    if pm.introspection.istype(varval, vartype, varsize)
        entry = varname + "=";
        delim = " ";
        vartype = lower(vartype);
        for i = 1 : varvalen
            if isa(varval(i), "cell")
                value = varval{i};
            else
                value = varval(i);
            end
            if strcmp(vartype, "string")
                entry = entry + pm.str.quote(string(value)) + delim;
            elseif strcmp(vartype, "integer")
                if rem(value, 1) == 0
                    entry = entry + value + delim;
                else
                    entry = entry + round(value) + delim;
                end
            elseif strcmp(vartype, "logical")
                if value
                    entry = entry + ".true." + delim;
                else
                    entry = entry + ".false." + delim;
                end
            elseif strcmp(vartype, "complex")
                entry = entry + "(" + string(real(value)) + "," + string(imag(value)) + ")" + delim;
            elseif strcmp(vartype, "real")
                entry = entry + value + delim;
            end
        end
    else
        disp(varname + " = ");
        disp(varval);
        error   ( newline ...
                + "The input " + varname + " specification value(s) displayed" + newline ...
                + "above must be conformable to a MATLAB " + vartype + " type," + newline ...
                + "with a maximum " + string(varsize) + " number of elements." + newline ...
                + "The specified value has the class:" + newline ...
                + newline ...
                + "    class(" + varname + ") = " + string(class(varval)) + newline ...
                + "with size:" + newline ...
                + newline ...
                + "    size(" + varname + ") = [" + join(string(size(varval)), ", ") + "]" + newline ...
                + newline ...
                );
    end
end
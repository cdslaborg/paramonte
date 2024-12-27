%>  \brief
%>  Return a Fortran-namelist-compatible conversion of the input ``varval``.<br>
%>
%>  \details
%>  This functionality is primarily used by the ParaMonte MATLAB internal
%>  routines to communicate information with Fortran shared libraries.<br>
%>  As such, it is of limited to most end users of the library.<br>
%>
%>  \param[in]  varname :   The input scalar MATLAB string containing the label to assign
%>                          to the namelist-converted value in the output ``entry``.<br>
%>                          The specified value of ``varname`` will be trimmed
%>                          (to remove leading and trailing blanks).<br>
%>  \param[in]  varval  :   The input value to be converted to namelist-compatible value.<br>
%>  \param[in]  vartype :   See the documentation of the corresponding
%>                          argument of [pm.introspection.istype()](@ref istype).<br>
%>  \param[in]  maxlen  :   See the documentation of the corresponding
%>                          argument of [pm.introspection.islenleq()](@ref islenleq).<br>
%>
%>  \return
%>  ``entry``           :   The output scalar MATLAB string containing the namelist-compatible
%>                          conversion of the input value ``varval`` and the given ``varname``
%>                          in the following format: ``varname=namelist-compatible-varval``.
%>
%>  \interface{getEntryNML}
%>  \code{.m}
%>
%>      entry = pm.fortran.getEntryNML(varname, varval, vartype, maxlen)
%>
%>  \endcode
%>
%>  \note
%>  If the input value is string, it will be quoted properly.<br>
%>  If the input ``varval`` is an array, its elements will be comma-separated.<br>
%>
%>  \see
%>  [pm.introspection.verify](@ref verify)<br>
%>  [pm.introspection.verified](@ref verified)<br>
%>  [pm.introspection.islenleq](@ref islenleq)<br>
%>  [pm.introspection.istype](@ref istype)<br>
%>
%>  \example{getEntryNML}
%>  \include{lineno} example/fortran/getEntryNML/main.m
%>  \output{getEntryNML}
%>  \include{lineno} example/fortran/getEntryNML/main.out.m
%>
%>  \final{getEntryNML}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:38 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function entry = getEntryNML(varname, varval, vartype, maxlen)
    varval = varval(:);
    varvalen = numel(varval);
    %varname = string(inputname(2));
    varname = string(strtrim(varname));
    if  pm.introspection.verified(varval, vartype, maxlen)
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
                entry = entry + pm.fortran.quote(string(value)) + delim;
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
            elseif strcmp(vartype, "real") || strcmpi(vartype, "float") || strcmpi(vartype, "single") || strcmpi(vartype, "double")
                entry = entry + value + delim;
            else
                help("pm.fortran.getEntryNML");
                disp("vartype");
                disp( vartype );
                error   ( newline ...
                        + "Unrecognized input ``vartype`` value." + newline ...
                        + newline ...
                        );
            end
        end
    else
        disp(varname + " = ");
        disp(varval);
        error   ( newline ...
                + "The input " + varname + " specification value(s) displayed" + newline ...
                + "above must be conformable to a MATLAB " + vartype + " type," + newline ...
                + "with a maximum " + string(maxlen) + " number of elements." + newline ...
                + "The specified value has the class:" + newline ...
                + newline ...
                + pm.io.tab() + "class(" + varname + ") = " + string(class(varval)) + newline ...
                + newline ...
                + "with size:" + newline ...
                + newline ...
                + pm.io.tab() + "size(" + varname + ") = [" + join(string(size(varval)), ", ") + "]" + newline ...
                + newline ...
                );
    end
end
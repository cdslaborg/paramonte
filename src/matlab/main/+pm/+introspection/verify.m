%>  \brief
%>  Verify the type and number of elements of the input ``varval``
%>  match the specified input ``vartype`` and ``varsize``.<br>
%>  Otherwise, halt the program by calling the MATLAB
%>  intrinsic ``error()`` with an appropriate
%>  descriptive error message.<br>
%>
%>  \details
%>  This function is merely a higher-level wrapper around
%>  the ParaMonte functionality [pm.introspection.istype()](@ref istype).<br>
%>  The goal is to offer a unified interface for argument verification and
%>  aborting the program with a suitable error message if
%>  the verification fails.<br>
%>
%>  \param[in]  varval  :   The input value to whose type and
%>                          number of elements is to be verified.
%>  \param[in]  vartype :   See the documentation of the corresponding
%>                          argument of [pm.introspection.istype()](@ref istype).
%>  \param[in]  varsize :   See the documentation of the corresponding
%>                          argument of [pm.introspection.istype()](@ref istype).
%>  \param[in]  varname :   The input scalar MATLAB string containing the label to assign
%>                          to the namelist-converted value in the output ``entry``.<br>
%>                          The specified value of ``varname`` will be trimmed
%>                          (to remove leading and trailing blanks) and used
%>                          only for outputting an error message.<br>
%>                          (**optional**, default = ``"value"``)
%>
%>  \interface{verify}
%>  \code{.m}
%>
%>      pm.introspection.verify(varval, vartype, varsize);
%>      pm.introspection.verify(varval, vartype, varsize, varname);
%>
%>  \endcode
%>
%>  \example{verify}
%>  \include{lineno} example/introspection/verify/main.m
%>  \output{verify}
%>  \include{lineno} example/introspection/verify/main.out.m
%>
%>  \final{verify}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:42 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function verify(varval, vartype, varsize, varname)
    if  nargin < 4
        varname = "value";
    else
        varname = string(strtrim(varname));
    end
    varval = varval(:);
    if ~pm.introspection.istype(varval, vartype, varsize)
        disp(varname + " = ");
        disp(varval);
        error   ( newline ...
                + "The input ``" + varname + "`` specification value(s) displayed" + newline ...
                + "above must be conformable to a MATLAB " + vartype + " type," + newline ...
                + "with a maximum " + string(varsize) + " number of elements." + newline ...
                + "The specified value has the class:" + newline ...
                + newline ...
                + pm.io.tab + "class(" + varname + ") = " + string(class(varval)) + newline ...
                + newline ...
                + "with size:" + newline ...
                + newline ...
                + pm.io.tab + "size(" + varname + ") = [" + join(string(size(varval)), ", ") + "]" + newline ...
                + newline ...
                );
    end
end
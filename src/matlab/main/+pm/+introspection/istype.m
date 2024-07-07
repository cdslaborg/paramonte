%>  \brief
%>  Return ``true`` if and only if the input ``varval`` conforms with the
%>  specified input type ``vartype`` and the specified maximum size ``varsize``.
%>  Otherwise, return ``false``.
%>
%>  \param[in]  varval  :   The input value to be converted to namelist-compatible value.
%>  \param[in]  vartype :   The input scalar MATLAB string containing the
%>                          expected type of the value given by the input ``varval``.<br>
%>                          The following type-conformance rules apply:<br>
%>                          <ol>
%>                              <li>    if ``vartype`` is ``"string"``, then ``varval`` can be
%>                                      either a MATLAB ``string`` or ``char``.<br>
%>                              <li>    if ``vartype`` is ``"integer"``, then ``varval`` can be
%>                                      either a MATLAB ``int8``, ``int16``, ``int32``, ``int64``,
%>                                      or a whole-number ``real`` value.<br>
%>                              <li>    if ``vartype`` is ``"logical"``, then ``varval`` can be
%>                                      either a MATLAB ``int8``, ``int16``, ``int32``, ``int64``,
%>                                      a MATLAB ``real``, or a MATLAB ``logical`` value.<br>
%>                              <li>    if ``vartype`` is ``"complex"``, then ``varval`` can be
%>                                      either a MATLAB ``int8``, ``int16``, ``int32``, ``int64``,
%>                                      a MATLAB ``real``, or a MATLAB ``complex`` value.<br>
%>                              <li>    if ``vartype`` is ``"real"``, then ``varval`` can be
%>                                      either a MATLAB ``int8``, ``int16``, ``int32``, ``int64``,
%>                                      or a MATLAB ``real`` value (e.g., ``float``, ``single``, ``double``).<br>
%>                                      For all other object types, the type-conformance is verified by
%>                                      passing the input ``varval`` and ``vartype`` directly to the
%>                                      MATLAB intrinsic function ``isa()``.<br>
%>                          </ol>
%>  \param[in]  varsize :   The input scalar MATLAB integer representing the
%>                          maximum allowed size of the input value ``varval``.<br>
%>                          (**optional**. If missing, the maximum length of the
%>                          input ``varval`` will not be checked.)
%>
%>  \return
%>  ``itis``            :   The output scalar MATLAB logical that is ``true`` if and only if
%>                          the input ``varval`` conforms with the specified input type ``vartype``
%>                          and the specified maximum size ``varsize``, otherwise, it is ``false``.<br>
%>
%>  \interface{istype}
%>  \code{.m}
%>
%>      itis = pm.introspection.istype(varval, vartype)
%>      itis = pm.introspection.istype(varval, vartype, varsize)
%>
%>  \endcode
%>
%>  \example{istype}
%>  \include{lineno} example/introspection/istype/main.m
%>  \output{istype}
%>  \include{lineno} example/introspection/istype/main.out.m
%>
%>  \final{istype}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:47 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = istype(varval, vartype, varsize)
    varvalen = numel(varval);
    itis = false;
    if  2 < nargin
        itis = varvalen <= varsize;
        if ~itis
            return;
        end
    end
    for i = 1 : varvalen
        if isa(varval(i), "cell")
            value = varval{i};
        else
            value = varval(i);
        end
        if  strcmpi(vartype, "string")
            itis = isa(value, "string") || isa(value, "char");
        elseif strcmpi(vartype, "integer")
            itis = isa(value, "int8") || isa(value, "int16") || isa(value, "int32") || isa(value, "int64");
            if ~itis && isreal(value)
                itis = rem(value, 1) == 0;
            end
        elseif strcmpi(vartype, "logical")
            itis = isa(value, "logical") || isreal(value) || isa(value, "int8") || isa(value, "int16") || isa(value, "int32") || isa(value, "int64");
        elseif strcmpi(vartype, "complex")
            itis = isnumeric(value);
        elseif strcmpi(vartype, "real") || strcmpi(vartype, "float") || strcmpi(vartype, "single") || strcmpi(vartype, "double")
            itis = isreal(value) || isa(value, "int8") || isa(value, "int16") || isa(value, "int32") || isa(value, "int64");
        else
            itis = isa(value, vartype);
        end
    end
end
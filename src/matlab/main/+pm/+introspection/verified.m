%>  \brief
%>  Return ``true`` if and only if the input ``varval`` conforms with the
%>  specified input type ``vartype`` and maximum length ``maxlen``,
%>  otherwise, return ``false``.<br>
%>
%>  \details
%>  If the input argument ``varval`` is a collection of values and the input ``vartype``
%>  is one of the five basic types (``string``, ``integer``, ``logical``, ``complex``, ``real``),
%>  then each element of the ``varval`` sequence will be tested against the input ``vartype``.<br>
%>
%>  **Beware** that this algorithm converts input values of type ``char`` to ``string``
%>  before further processing and type checking, that is, the types ``char`` to ``string``
%>  are assumed to be conformable and compatible with each other, like most other sane languages.<br>
%>
%>  \param[in]  varval  :   The input value whose type and length must be verified.
%>  \param[in]  vartype :   The input scalar MATLAB string containing the
%>                          expected type of the value given by the input ``varval``.<br>
%>                          The following type-conformance rules apply:<br>
%>                          <ol>
%>                              <li>    if ``vartype`` is ``"string"``, then ``varval`` can be
%>                                      either a MATLAB ``string`` or ``char``.<br>
%>                                      An input value of type ``char`` is always
%>                                      converted to ``string`` before further processing.<br>
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
%>  @param[in]  maxlen  :   The input positive scalar MATLAB integer representing the
%>                          maximum allowed size of the input value ``varval``.<br>
%>                          The number of elements of ``varval`` is measured by the
%>                          function [pm.introspection.islenleq()](@ref islenleq)<br>
%>
%>  \return
%>  ``itis``            :   The output scalar MATLAB logical that is ``true`` if and only if
%>                          the input ``varval`` conforms with the specified input type ``vartype``
%>                          and maximum length, otherwise, it is ``false``.<br>
%>
%>  \interface{verified}
%>  \code{.m}
%>
%>      itis = pm.introspection.verified(varval, vartype, maxlen)
%>
%>  \endcode
%>
%>  \see
%>  [pm.introspection.verify](@ref verify)<br>
%>  [pm.introspection.verified](@ref verified)<br>
%>  [pm.introspection.islenleq](@ref islenleq)<br>
%>  [pm.introspection.istype](@ref istype)<br>
%>
%>  \example{verified}
%>  \include{lineno} example/introspection/verified/main.m
%>  \output{verified}
%>  \include{lineno} example/introspection/verified/main.out.m
%>
%>  \final{verified}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:47 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = verified(varval, vartype, maxlen)
    itis = pm.introspection.islenleq(varval, maxlen) && pm.introspection.istype(varval, vartype);
end
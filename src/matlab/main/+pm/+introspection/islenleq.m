%>  \brief
%>  Return ``true`` if and only if the input ``varval``
%>  has the specified maximum size ``maxlen``,
%>  Otherwise, return ``false``.
%>
%>  \details
%>  This routine is primarily used within the procedures of
%>  the ParaMonte library for type introspection and verification
%>  in high-level programming languages.<br>
%>
%>  @param[in]  varval  :   The input value whose length is to verified to
%>                          be less than or equal to the input value ``maxlen``.<br>
%>                          Beware that,
%>                          <ol>
%>                              <li>    an input scalar string has a length of ``1``.<br>
%>                              <li>    an input scalar value has a length of ``1``.<br>
%>                              <li>    the length of a non-scalar object is measured
%>                                      by the MATLAB intrinsic ``numel()``.<br>
%>                          </ol>
%>  @param[in]  maxlen  :   The input positive scalar MATLAB integer representing the
%>                          maximum allowed size of the input value ``varval``.<br>
%>
%>  @return
%>  ``itis``            :   The output scalar MATLAB logical that is ``true`` if and only if
%>                          the input ``varval`` has a length that is less than or equal to
%>                          the specified maximum size ``maxlen``, otherwise, it is ``false``.<br>
%>
%>  \interface{islenleq}
%>  \code{.m}
%>
%>      itis = pm.introspection.islenleq(varval, maxlen)
%>
%>  \endcode
%>
%>  \see
%>  [pm.introspection.verify](@ref verify)<br>
%>  [pm.introspection.verified](@ref verified)<br>
%>  [pm.introspection.islenleq](@ref islenleq)<br>
%>  [pm.introspection.istype](@ref istype)<br>
%>
%>  \example{islenleq}
%>  \include{lineno} example/introspection/islenleq/main.m
%>  \output{islenleq}
%>  \include{lineno} example/introspection/islenleq/main.out.m
%>
%>  \final{islenleq}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:47 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = islenleq(varval, maxlen)
    if  ischar(varval)
        varval = string(varval);
    end
    itis = numel(varval) <= maxlen;
end
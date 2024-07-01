%>  \brief
%>  Return a scalar MATLAB whole-number containing
%>  the length of the input scalar or vector object.
%>
%>  \details
%>  The following rules apply to computing the length:
%>  <ol>
%>      <li>    If the input ``obj`` is a scalar object of type ``char``
%>              or ``string`` and empty, then the output ``val`` is ``0``.
%>
%>      <li>    If the input ``obj`` is a scalar object of type ``char``
%>              or ``string`` and non-empty, then the output ``val`` is ``1``.
%>
%>      <li>    If the input ``obj`` is a scalar object of type other
%>              than ``char`` or ``string``, then the output ``val`` is ``1``.
%>              This is similar to the behavior of the MATLAB intrinsic ``length()``.
%>
%>      <li>    If the input ``obj`` is a vector and empty, then the output ``val`` is ``0``.
%>              This is similar to the behavior of the MATLAB intrinsic ``length()``.
%>
%>      <li>    If the input ``obj`` is a vector and non-empty,
%>              then the output ``val`` is the vector length as returned
%>              by ``length(obj)`` minus the number of empty elements.
%>  </ol>
%>
%>  \param[in]  obj :   The input scalar or vector object whose length
%>                      is to returned according to the rules above.
%>
%>  \return
%>   `val`          :   The output scalar MATLAB whole-number
%>                      representing the length of the input object.
%>
%>  \interface{len}
%>  \code{.m}
%>
%>      val = pm.array.len(obj);
%>
%>  \endcode
%>
%>  \example{len}
%>  \include{lineno} example/array/len/main.m
%>  \output{len}
%>  \include{lineno} example/array/len/main.out.m
%>
%>  \final{len}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:25 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function val = len(obj)
    if isa(obj, "char")
        obj = string(obj);
    end
    % First ensure the input ``obj`` is a vector.
    if numel(obj) ~= length(obj)
        error   ( newline ...
                + "The input object must be a vector." + newline ...
                + newline ...
                );
    end
    val = 0;
    for i = 1 : length(obj)
        try
            if ~strcmp(string(obj(i)), "")
                val = val + 1;
            end
        catch
            val = val + 1;
        end
    end
end
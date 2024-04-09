function val = len(obj)
    %
    %   Return a scalar MATLAB whole-number containing
    %   the length of the input scalar or vector object.
    %
    %   The following rules apply to computing the length:
    %
    %       -   If the input ``obj`` is a scalar object of type ``char``
    %           or ``string`` and empty, then the output ``val`` is ``0``.
    %
    %       -   If the input ``obj`` is a scalar object of type ``char``
    %           or ``string`` and non-empty, then the output ``val`` is ``1``.
    %
    %       -   If the input ``obj`` is a scalar object of type other
    %           than ``char`` or ``string``, then the output ``val`` is ``1``.
    %           This is similar to the behavior of the MATLAB intrinsic ``length()``.
    %
    %       -   If the input ``obj`` is a vector and empty, then the output ``val`` is ``0``.
    %           This is similar to the behavior of the MATLAB intrinsic ``length()``.
    %
    %       -   If the input ``obj`` is a vector and non-empty,
    %           then the output ``val`` is the vector length as returned
    %           by ``length(obj)`` minus the number of empty elements.
    %
    %   Parameters
    %   ----------
    %
    %       obj
    %
    %           The input scalar or vector object whose length
    %           is to returned according to the rules above.
    %
    %   Returns
    %   -------
    %
    %       val
    %
    %           The output scalar MATLAB whole-number
    %           representing the length of the input object.
    %   Interface
    %   ---------
    %
    %       pm.array.len(obj)
    %
    %   Example
    %   -------
    %
    %       pm.array.len(1) % = 1
    %       pm.array.len("") % = 0
    %       pm.array.len([]) % = 0
    %       pm.array.len('paramonte') % = 1
    %       pm.array.len("paramonte") % = 1
    %       pm.array.len(["pm", 'array']) % = 2
    %       pm.array.len(["pm", 'array', []]) % = 2
    %
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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
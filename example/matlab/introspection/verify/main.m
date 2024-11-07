cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.introspection.verify(10000, "integer", 1, "varname");
pm.introspection.verify(10000, "integer", 1);
pm.introspection.verify("text", "string", 1);
pm.introspection.verify('text', "char", 4);
pm.introspection.verify('text', "char", 5);
pm.introspection.verify(1, "complex", 1);
pm.introspection.verify(1, "real", 1);
pm.introspection.verify(int32(1), "int32", 1);

try
    pm.introspection.verify(10000, "string", 1);
catch me
    warning(string(me.identifier) + string(me.message));
end
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('from = struct("key1", 1, "key2", "val2", "Key2", "val2duplicate", "key3", false, "key4", []);')
                from = struct("key1", 1, "key2", "val2", "Key2", "val2duplicate", "key3", false, "key4", []);
pm.matlab.show('to = struct("key1", "val1", "Key2", "Val2");')
                to = struct("key1", "val1", "Key2", "Val2");

pm.matlab.show()
pm.matlab.show('new = pm.matlab.copy(from);')
                new = pm.matlab.copy(from);
pm.matlab.show('new')
pm.matlab.show( new )

pm.matlab.show()
pm.matlab.show('new = pm.matlab.copy(from, [], "key3");')
                new = pm.matlab.copy(from, [], "key3");
pm.matlab.show('new')
pm.matlab.show( new )

pm.matlab.show()
pm.matlab.show('new = pm.matlab.copy(from, to, "key3");')
                new = pm.matlab.copy(from, to, "key3");
pm.matlab.show('new')
pm.matlab.show( new )

pm.matlab.show()
pm.matlab.show('new = pm.matlab.copy(from, to);')
                new = pm.matlab.copy(from, to);
pm.matlab.show('new')
pm.matlab.show( new )
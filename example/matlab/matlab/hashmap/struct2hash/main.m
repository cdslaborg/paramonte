cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('s = struct("key1", 1, "key2", "val2", "Key2", "val2duplicate", "key3", false, "key4", []);')
                s = struct("key1", 1, "key2", "val2", "Key2", "val2duplicate", "key3", false, "key4", []);

pm.matlab.show()
pm.matlab.show('hashmap = pm.matlab.hashmap.struct2hash(s);')
                hashmap = pm.matlab.hashmap.struct2hash(s);
pm.matlab.show('hashmap')
pm.matlab.show( hashmap )

pm.matlab.show()
pm.matlab.show('hashmap = pm.matlab.hashmap.struct2hash(s, "key1");')
                hashmap = pm.matlab.hashmap.struct2hash(s, "key1");
pm.matlab.show('hashmap')
pm.matlab.show( hashmap )

pm.matlab.show()
pm.matlab.show('hashmap = pm.matlab.hashmap.struct2hash(s, [], true);')
                hashmap = pm.matlab.hashmap.struct2hash(s, [], true);
pm.matlab.show('hashmap')
pm.matlab.show( hashmap )

pm.matlab.show()
pm.matlab.show('hashmap = pm.matlab.hashmap.struct2hash(s, "key1", true);')
                hashmap = pm.matlab.hashmap.struct2hash(s, "key1", true);
pm.matlab.show('hashmap')
pm.matlab.show( hashmap )

pm.matlab.show()
pm.matlab.show('hashmap = pm.matlab.hashmap.struct2hash(s, "key1", [], true);')
                hashmap = pm.matlab.hashmap.struct2hash(s, "key1", [], true);
pm.matlab.show('hashmap')
pm.matlab.show( hashmap )
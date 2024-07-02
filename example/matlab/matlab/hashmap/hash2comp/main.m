cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('hashmap = {"key1", 2, "key2", "hash", "key3", true, "key4", 0};')
                hashmap = {"key1", 2, "key2", "hash", "key3", true, "key4", 0};
pm.matlab.show('object = struct("key1", 1, "key2", "val2", "Key2", "val2duplicate", "key3", false, "key4", []);')
                object = struct("key1", 1, "key2", "val2", "Key2", "val2duplicate", "key3", false, "key4", []);

pm.matlab.show()
pm.matlab.show('objnew = pm.matlab.hashmap.hash2comp(hashmap, object);')
                objnew = pm.matlab.hashmap.hash2comp(hashmap, object);
pm.matlab.show('objnew')
pm.matlab.show( objnew )
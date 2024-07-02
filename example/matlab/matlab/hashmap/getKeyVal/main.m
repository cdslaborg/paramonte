cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('hashmap = {"key1", 1, "key2", "val2", "key3", false};')
                hashmap = {"key1", 1, "key2", "val2", "key3", false};

pm.matlab.show()
pm.matlab.show('[val, failed] = pm.matlab.hashmap.getKeyVal("key2", hashmap);')
                [val, failed] = pm.matlab.hashmap.getKeyVal("key2", hashmap);
pm.matlab.show('val')
pm.matlab.show( val )
pm.matlab.show('failed')
pm.matlab.show( failed )
assert(val == "val2" && ~failed)

pm.matlab.show()
pm.matlab.show('[val, failed] = pm.matlab.hashmap.getKeyVal("key3", hashmap);')
                [val, failed] = pm.matlab.hashmap.getKeyVal("key3", hashmap);
pm.matlab.show('val')
pm.matlab.show( val )
pm.matlab.show('failed')
pm.matlab.show( failed )
assert(val == false && ~failed)

pm.matlab.show()
pm.matlab.show('[val, failed] = pm.matlab.hashmap.getKeyVal("key3", hashmap(1:4));')
                [val, failed] = pm.matlab.hashmap.getKeyVal("key3", hashmap(1:4));
pm.matlab.show('val')
pm.matlab.show( val )
pm.matlab.show('failed')
pm.matlab.show( failed )
assert(isempty(val) && failed)

pm.matlab.show()
pm.matlab.show('[val, failed] = pm.matlab.hashmap.getKeyVal("key2", hashmap(1:4));')
                [val, failed] = pm.matlab.hashmap.getKeyVal("key2", hashmap(1:4));
pm.matlab.show('val')
pm.matlab.show( val )
pm.matlab.show('failed')
pm.matlab.show( failed )
assert(val == "val2" && ~failed)

pm.matlab.show()
pm.matlab.show('try; pm.matlab.hashmap.getKeyVal("key2", hashmap(1:3)); catch me; disp(string(me.identifier) + string(me.message)); end')
                try; pm.matlab.hashmap.getKeyVal("key2", hashmap(1:3)); catch me; disp(string(me.identifier) + string(me.message)); end
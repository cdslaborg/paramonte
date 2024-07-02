cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('hashmap = {"key1", 1, "key2", "val2", "key3", false, "key2", "last"};')
                hashmap = {"key1", 1, "key2", "val2", "key3", false, "key2", "last"};

pm.matlab.show()
pm.matlab.show('[keyval, hashout] = pm.matlab.hashmap.popKeyVal("key2", hashmap);')
                [keyval, hashout] = pm.matlab.hashmap.popKeyVal("key2", hashmap);
pm.matlab.show('keyval')
pm.matlab.show( keyval )
pm.matlab.show('hashout')
pm.matlab.show( hashout )

pm.matlab.show()
pm.matlab.show('[keyval, hashout] = pm.matlab.hashmap.popKeyVal("key3", hashmap);')
                [keyval, hashout] = pm.matlab.hashmap.popKeyVal("key3", hashmap);
pm.matlab.show('keyval')
pm.matlab.show( keyval )
pm.matlab.show('hashout')
pm.matlab.show( hashout )

pm.matlab.show()
pm.matlab.show('[keyval, hashout] = pm.matlab.hashmap.popKeyVal("key3", hashmap(1:4));')
                [keyval, hashout] = pm.matlab.hashmap.popKeyVal("key3", hashmap(1:4));
pm.matlab.show('keyval')
pm.matlab.show( keyval )
pm.matlab.show('hashout')
pm.matlab.show( hashout )

pm.matlab.show()
pm.matlab.show('[keyval, hashout] = pm.matlab.hashmap.popKeyVal({"key2", "key1"}, hashmap(1:4));')
                [keyval, hashout] = pm.matlab.hashmap.popKeyVal({"key2", "key1"}, hashmap(1:4));
pm.matlab.show('keyval')
pm.matlab.show( keyval )
pm.matlab.show('hashout')
pm.matlab.show( hashout )

pm.matlab.show()
pm.matlab.show('try; [keyval, hashout] = pm.matlab.hashmap.popKeyVal("key2", hashmap(1:3)); catch me; disp(string(me.identifier) + string(me.message)); end')
                try; [keyval, hashout] = pm.matlab.hashmap.popKeyVal("key2", hashmap(1:3)); catch me; disp(string(me.identifier) + string(me.message)); end
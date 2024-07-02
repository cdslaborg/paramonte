cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('hashmap = {"key1", 1, "key2", "val2"};')
                hashmap = {"key1", 1, "key2", "val2"};

pm.matlab.show()
pm.matlab.show('hashnew = pm.matlab.hashmap.addKeyVal("key3", false, {});')
                hashnew = pm.matlab.hashmap.addKeyVal("key3", false, {});
pm.matlab.show('hashnew')
pm.matlab.show( hashnew )

pm.matlab.show()
pm.matlab.show('hashnew = pm.matlab.hashmap.addKeyVal("key2", "val2", hashmap);')
                hashnew = pm.matlab.hashmap.addKeyVal("key2", "val2", hashmap);
pm.matlab.show('hashnew')
pm.matlab.show( hashnew )
assert(all([hashnew{:}] == [hashmap{:}]))
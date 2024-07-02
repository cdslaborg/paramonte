cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('hashmap = {"key1", 1, "key2", "val2"};')
                hashmap = {"key1", 1, "key2", "val2"};

pm.matlab.show()
pm.matlab.show('hashnew = pm.matlab.hashmap.repKeyVal("key2", 2, hashmap);')
                hashnew = pm.matlab.hashmap.repKeyVal("key2", 2, hashmap);
pm.matlab.show('hashnew')
pm.matlab.show( hashnew )

pm.matlab.show()
pm.matlab.show('hashnew = pm.matlab.hashmap.repKeyVal("key3", true, hashmap);')
                hashnew = pm.matlab.hashmap.repKeyVal("key3", true, hashmap);
pm.matlab.show('hashnew')
pm.matlab.show( hashnew )

pm.matlab.show()
pm.matlab.show('hashnew = pm.matlab.hashmap.repKeyVal("key3", false, hashmap);')
                hashnew = pm.matlab.hashmap.repKeyVal("key3", false, hashmap);
pm.matlab.show('hashnew')
pm.matlab.show( hashnew )
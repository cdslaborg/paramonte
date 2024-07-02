cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show('weblinks = pm.lib.weblinks()')
                weblinks = pm.lib.weblinks()
pm.matlab.show()

pm.matlab.show()
pm.matlab.show('pm.matlab.show(weblinks)')
                pm.matlab.show(weblinks)

pm.matlab.show()
pm.matlab.show('pm.matlab.show(weblinks, [], true)')
                pm.matlab.show(weblinks, [], true)

pm.matlab.show()
pm.matlab.show('pm.matlab.show(pm.lib.weblinks(), "weblinks", true)')
                pm.matlab.show(pm.lib.weblinks(), "weblinks", true)

pm.matlab.show()
pm.matlab.show('pm.matlab.show(pm.lib.weblinks(), "weblinks")')
                pm.matlab.show(pm.lib.weblinks(), "weblinks")

pm.matlab.show()
pm.matlab.show('struct("key1", "val1", "Key2", "Val2"), "s"')
pm.matlab.show( struct("key1", "val1", "Key2", "Val2"), "s" )

pm.matlab.show()
pm.matlab.show('{"key1", 1, "key2", "val2"}, "mycell"')
pm.matlab.show( {"key1", 1, "key2", "val2"}, "mycell" )

pm.matlab.show()
pm.matlab.show('[1, 2, 3, 4]', "vec")
pm.matlab.show( [1, 2, 3, 4] , "vec")

pm.matlab.show()
pm.matlab.show('"string", "str"')
pm.matlab.show( "string", "str" )
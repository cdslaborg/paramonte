cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('cd(tempdir); pwd')
                cd(tempdir); pwd

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("File.Ext")')
pm.matlab.show( pm.sys.path.abs("File.Ext") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("..\File.Ext")')
pm.matlab.show( pm.sys.path.abs("..\File.Ext") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("..\..\File.Ext")')
pm.matlab.show( pm.sys.path.abs("..\..\File.Ext") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs(".\File.Ext")')
pm.matlab.show( pm.sys.path.abs(".\File.Ext") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("*.txt")')
pm.matlab.show( pm.sys.path.abs("*.txt") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("..")')
pm.matlab.show( pm.sys.path.abs("..") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("..\..\..")')
pm.matlab.show( pm.sys.path.abs("..\..\..") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("Folder\")')
pm.matlab.show( pm.sys.path.abs("Folder\") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("D:\A\..\B")')
pm.matlab.show( pm.sys.path.abs("D:\A\..\B") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs("\\Server\Folder\Sub\..\File.ext")')
pm.matlab.show( pm.sys.path.abs("\\Server\Folder\Sub\..\File.ext") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs(["..", "new"])')
pm.matlab.show( pm.sys.path.abs(["..", "new"]) )

pm.matlab.show()
pm.matlab.show("pm.sys.path.abs({'..', 'new'})")
pm.matlab.show( pm.sys.path.abs({'..', 'new'}) )

pm.matlab.show()
pm.matlab.show('pm.sys.path.abs(".", "fat")')
pm.matlab.show( pm.sys.path.abs(".", "fat") )
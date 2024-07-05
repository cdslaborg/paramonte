cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob("https://apod.nasa.gov/apod/astropix.html")')
pm.matlab.show( pm.sys.path.webglob("https://apod.nasa.gov/apod/astropix.html") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob("https://github.com/cdslaborg/paramontex/blob/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/sampling_multivariate_normal_distribution_via_paradram/out/mvn_serial_process_1_report.txt")')
pm.matlab.show( pm.sys.path.webglob("https://github.com/cdslaborg/paramontex/blob/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/sampling_multivariate_normal_distribution_via_paradram/out/mvn_serial_process_1_report.txt") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob("*.m")')
pm.matlab.show( pm.sys.path.webglob("*.m") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob("example/array/**")')
pm.matlab.show( pm.sys.path.webglob("example/array/**") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob("example/array/**.m")')
pm.matlab.show( pm.sys.path.webglob("example/array/**.m") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob("example/array/*.m")')
pm.matlab.show( pm.sys.path.webglob("example/array/*.m") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.webglob(".*")')
pm.matlab.show( pm.sys.path.webglob(".*") )

pm.matlab.show()
pm.matlab.show('[paths, isdir] = pm.sys.path.webglob("**"); paths(~isdir) % get all files in directory tree.')
                [paths, isdir] = pm.sys.path.webglob("**"); paths(~isdir) % get all files in directory tree.

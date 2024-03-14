REM Add the library path to the environment PATH variable.
setlocal EnableDelayedExpansion
set "PATH=..\..\..\lib;%PATH%"
REM Generate build, build, and run
cmake . -G "NMake Makefiles" && nmake && nmake run
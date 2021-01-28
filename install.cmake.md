The easiest way to build and install the ParaMonte library and its examples is via the Linux/macOS 
Bash and Windows Batch install scripts that exist at the root of the ParaMonte repository on GitHub. 
The scripts will automatically install the library and run the tests (if requested) and the examples provided.  

However, if you still really prefer to build library directly via the cmake software and 
manually go through all the rest of the build, testing, and installation processes, 
you must follow the guidelines below to predefine the various aspects of the library and its tests.  

**Note** that the ParaMonte library build via cmake is exclusively available on Linux, macOS, and Microsoft WSL platforms.  

1.  Ensure a [ParaMonte-compatible](CHANGES.md) version of the GNU or Intel compiler suite 
    is already installed on our system and the paths to different compiler components and wrappers 
    already exist in the environmental `PATH` variable. 
    (**This step would be automatically done if you used the `install.sh` script.**)  
1.  Ensure a ParaMonte-compatible cmake installation (`3.14.0` or newer) exists on your system 
    and the path to the cmake executable already exist in the environmental `PATH` variable. 
    (**This step would be automatically done if you used the `install.sh` script.**)  
1.  Navigate to the root directory of the library (where the `.git` file exists).  
1.  Then, assuming that you are using the GNU compiler suite, you can configure the library's build minimally via,  
    ```bash  
    mkdir build && \
    cd build && \
    cmake -DPMCS=GNU ..
    ```  
    (**This step would be automatically done if you used the `install.sh` script.**)  
1.  The above command will configure the library with the default build settings. 
    To change the default settings, the following macros can be also defined,  
    ```bash  
    mkdir build && \
    cd build && \
    cmake \
    --verbose=1 \
    -DPMCS="the compiler suite: GNU/INTEL" \
    -DPMLIB_NAME="you_must_follow_the_guidelines_here:https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds" \
    -DCMAKE_Fortran_COMPILER="/path/to/the/specific/fortran/compiler/executable/to/be/used/example:/usr/bin/gfortran-10" \
    -DMATLAB_ROOT_DIR="/if/building/for/matlab/provide/path/to/the/matlab/root/installation/directory/where/bin/dir/exists" \
    -DMPIEXEC_EXECUTABLE="/path/to/the/specific/mpi/launcher/to/be/used" \
    -DINTERFACE_LANGUAGE="the single programming language interface: c/cpp/fortran/matlab/python/... (default = c)" \
    -DMPI_ENABLED="whether building for MPI parallelization: true/false (default = false)" \
    -DCAFTYPE="The Coarray parallelization type: none/single/shared/distributed (default = none)" \
    -DBTYPE="library build: debug/testing/release (default = release)" \
    -DLTYPE="library type: static/shared (default = static)" \
    -DHEAP_ARRAY_ENABLED="whether to perform all memory allocations on the heap: true/false (default = false)" \
    -DCFI_ENABLED="whether C-Fortran interoperation must be enabled: true/false (default = true)" \
    -DOMP_ENABLED="whether OpenMP parallelization must be enabled: true/false (default = false)" \
    -DMPILIB_NAME="The MPI library name: impi/mpich/openmpi (impi stands for Intel MPI)" \
    -DSAMPLER_TEST_ENABLED="whether the library's sampler tests should be enabled: true/false (default = false)" \
    -DBASIC_TEST_ENABLED="whether the library's basic tests should be enabled: true/false (default = false)" \
    -DPERFPROF_ENABLED="whether the library should be build for performance profiling: true/false (default = false)" \
    -DCODECOV_ENABLED="whether the library should be build for code coverage: true/false (default = false)" \
    -DOS_IS_WSL="whether the current operating system is Microsoft WSL" \
    ..
    ```  
    (**All of the above configuration settings could be easily done if you used the `install.sh` script.**)  
1.  Once you configure the library's build, try,  
    ```bash  
    make -j 3 && make install
    ```  
    Replace the number `3` with as many number of cores you wish to use to build the library. 
    (**This step would be automatically done if you used the `install.sh` script.**)  
1.  The ParaMonte library should be now built and installed in the `lib` directory. 
    If you built the library for testing or code coverage purposes, then a `testParaMonte` 
    executable must also exist in the subdirectory `test/bin/`. The onus is on you, the developer/user, to ensure 
    all runtime libraries and dependencies can be found.  
    (**All dependency checks and test runs would be automatically done if you used the `install.sh` script.**)  
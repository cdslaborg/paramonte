# ParaMonte release notes  

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases).  

>   To learn about language-specific changes, visit:  
>
>   +   [ParaMonte C CHANGES.md](./src/c/CHANGES.md).  
>   +   [ParaMonte C++ CHANGES.md](./src/cpp/CHANGES.md).  
>   +   [ParaMonte Fortran CHANGES.md](./src/fortran/CHANGES.md).  
>   +   [ParaMonte MATLAB CHANGES.md](./src/matlab/CHANGES.md).  

## **Version 2.x.x**  

### **Version 2.0.0** -- March 8, 2024 (pre-release)  

**Major Enhancements**  

+   This is a major interface-breaking ParaMonte library release:  
    +   The ParaMonte C and C++ have considerably changed with new interfaces and functionalities. 
    +   The ParaMonte Fortran library has grown extensively into a full-blown Machine Learning library.  
    +   The ParaMonte MATLAB and Python libraries have also grown extensively with modernized simple interfaces.  

The current pre-release contains only the ParaMonte C, C++, Fortran, and MATLAB implementations.  
Work and significant testing is underway to release the ParaMonte Python and R implementations.  

**Compiler Compatibility**  

| Compiler Suite                    | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| GNU Compiler Collection >= 10.3   | ✅              | ✅            | ✅             | ✅            |  
| Intel OneAPI >= 2021.8.0          | ✅              | ✅            | ✅             | ❌            |  

**Essential Dependency Compatibility**  

| Dependency                        | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| CMake >= 3.21                     | ✅              | ✅            | ✅             | ✅            |  

**Optional Dependency Compatibility**  

| Dependency                        | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| Intel MPI (IMPI) >= 2021.8        | ✅              | ✅            | ✅             | ❌            |  
| MPICH MPI (MMPI) >= 3             | ❌              | ✅            | ✅             | ✅            |  
| OpenMPI (OMPI) >= 4               | ❌              | ✅            | ✅             | ✅            |  

## **Version 1.x.x**  

### **Version 1.5.1** -- January 1, 2021  

**Minor Enhancements**  

+   This is a minor enhancement release but is a major step 
    toward further portability of the ParaMonte library. All ParaMonte 
    library dependencies are now properly handled and recognized at runtime
    without such aggressive actions as permanently redefining the environmental
    path variables, most importantly, `PATH` and `LD_LIBRARY_PATH` on Linux/macOS.

+   The ParaMonte routines are now capable of handling user-input 
    file paths that contain white space (blank) or other exotic characters.

+   The Bash build script for the ParaMonte C/C++/Fortran examples can now handle
    file paths that contain white-space (blank) or other exotic characters.

+   The shared (dynamic) library file naming convention is now changed from
    using "dynamic" in the file name to using "shared" in the file name.

+   Improved error-handling messages in the examples' build and run scripts.

+   Typo-fixes in the documentation of the library.

**Compiler Compatibility**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 8.4     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 19.1.1    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.1.1.216 Build 20200306 / Intel(R) MPI Library 2019 Update 7 for Windows  
+   **Linux**: Intel Parallel Studio Version 19.1.1.217 20200306 / Intel(R) MPI Library for Linux OS, Version 2019 Update 7 Build 20200312
+   **Linux**: GNU 10.1.0 / Open-MPI 4.0.3  
+   **Linux**: GNU 10.1.0 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  
+   **macOS**: GNU 10.2.0 / MPICH 3.3.2  

### **Version 1.5.0** -- December 17, 2020

**Major Enhancements**  

+   This version introduces numerous performance and accuracy enhancements to the ParaMonte library.
+   The entire ParaMonte library is now fully documented and verified with over 866 tests that cover close 
    to 100% of all lines and functions in the library.
+   New prebuilt libraries with GNU compilers and Open-MPI on Linux are added.
+   New flags are now added to the build scripts of the library that automate the process of code coverage generation.
+   The `testing` builds are now removed from the ParaMonte release page as this build is mostly useful for development purposes.  
+   The issue of Windows file locking, that led to the occasional crashes of the 
    ParaDRAM and ParaDISE simulations in `multiChain` parallelism mode, is now resolved.

**Minor Enhancements**  

+   Minor enhancements to the ParaMonte C/C++/Fortran example build scripts `build.sh` and `build.bat`.  
+   The default build settings are now limited to `heap` memory allocation with `dynamic` library builds
    for only serial and MPI parallelization for all languages. The `stack` memory allocation results in 
    a ~10% gain in the efficiency of the code. The benefits of stack-memory builds are marginal and are 
    often problematic, in particular, for usage with non-compiled languages. However, users can still 
    build the library with `stack` memory allocation by specifying the appropriate build flags with
    the `install.sh` on Unix or `install.bat` script on Windows systems. For further information, 
    see the installation guidelines on the ParaMonte documentation website.
+   All temporary array creations in debug mode are now resolved, 
    except when Intel compilers are used, in which case, the debug warning messages are silenced.

**Compiler Compatibility**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 8.4     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 19.1.1    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.1.1.216 Build 20200306 / Intel(R) MPI Library 2019 Update 7 for Windows  
+   **Linux**: Intel Parallel Studio Version 19.1.1.217 20200306 / Intel(R) MPI Library for Linux OS, Version 2019 Update 7 Build 20200312
+   **Linux**: GNU 10.1.0 / Open-MPI 4.0.3  
+   **Linux**: GNU 10.1.0 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  
+   **macOS**: GNU 10.2.0 / MPICH 3.3.2  

### **Version 1.4.1** -- November 15, 2020

**Enhancements**  

+   The ParaMonte C/C++/Fortran example build scripts `build.sh` and `build.bat` 
    now check for the existence of both Intel and GNU compilers in the appropriate 
    order that is automatically inferred at the compilation time. Also, the dependencies 
    on the MPI compiler wrappers is now removed as the MPI libraries are not required to 
    build the ParaMonte examples, even in cases of parallel ParaMonte example builds.

**Compiler support**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 10.1.0 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.4.0** -- October 29, 2020

**Enhancements**  

+   The IO debugging info of all ParaMonte samplers have been enhanced. 
    In cases of wrong syntax or syntax-breaking input values in the simulation 
    output files, the error messages are now more informative and point directly 
    to the exact location of error in the input file.  

+   The Integrated Autocorrelation (ACT) for sample refinement in ParaDRAM 
    sampler of ParaMonte is now set to the average of all variables' ACT values 
    instead of the maximum ACT value. This will lead to less aggressive decorrelation 
    of the final sample, which means significantly larger final sample sizes, without 
    compromising the i.i.d. property of the final refined sample. This behavior can 
    be reversed back to the original by specifying "max" or "maximum" along with 
    the requested refinement method, `SampleRefinementMethod = "batchmeans max"` 
    or `SampleRefinementMethod = "BatchMeans-max"` (case-insensitive).

**Compiler support**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 9.1 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.3.0** -- October 3, 2020

**Enhancements**  

+   A new simulation specification `overwriteRequested` has 
    been added to all ParaMonte samplers. If TRUE and the 
    ParaMonte sampler detects an existing set of old simulation 
    output files in the output path of the current simulation with 
    the same names as the output file names of the current simulation, 
    then, the ParaMonte sampler will overwrite the existing simulation files.  

**Compiler support**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 9.1 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.2.0** -- September 22, 2020

**Enhancements**  

+   The post-processing report in the output report file 
    of ParaDRAM simulation has been significantly improved:  
    +   The parallel simulation summary now also provides the 
        predicted strong-scaling speedup behavior of the parallel 
        ParaDRAM simulations in "single chain" parallelism mode. 
        This can help make wiser decisions regarding the number 
        of processors for similar parallel simulations in the future.  

+   The ParaDRAM restart output file in ASCII mode now contains all 
    proposal updates, including the first user-specified proposal specs.

+   All parallel simulations now avoid the unnecessary creation of 
    temporary files by all processors for System and OS operations. 
    This is particularly important for large-scale parallel simulations.
    As a side effect, this will also potentially improve the runtime 
    performances of the simulation.  

+   Major enhancements has been made to the parallel simulation 
    performance analysis reported in the post-processing section 
    of the ParaDRAM simulation output `_report.txt` files.  

+   The ParaMonte C/C++ build processes are now separate from each other.  

**Bug fixes**  

+   minor typo fixes

**Compiler support**  

| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 9.1 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.1.0** -- May 27, 2020  

**New features**  

+   The MATLAB interface to the ParaMonte library is now ready to use, in addition 
    to the existing C/C++/Fortran/Python Programming language interfaces to ParaMonte.  

**Enhancements**  

+   The simulation summary report in the output report file has been improved.
+   The error handling in the build Batch-scripts on Windows is now greatly improved.

**Bug fixes**  

+   The build Batch-script for the ParaMonte examples on Windows now properly builds and runs coarray applications in parallel.
+   The fully-deterministic restart functionality is now functional also when outputChainFormat="verbose" in ParaDRAM simulations.

**Compiler support**  

| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 8.3 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.166 20191121  
+   **macOS**: GNU 9.3 / Open-MPI 4.0  

### **Version 1.0.0** -- January 1, 2020 -- Initial release  

This is the first public release of the ParaMonte library.  

**New features**  

+   ParaDRAM sampler: **Para**llel **D**elayed-**R**ejection **A**daptive Metropolis-Hastings **M**arkov Chain Monte Carlo Sampler.  
+   ParaMonte Interface to the C/C++/Fortran/Python Programming languages.  
+   ParaMonte simulation-output visualization via the ParaMonte Python interface.  

**Compiler support**  

| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 8.3 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.166 20191121  
+   **macOS**: GNU 9.3 / Open-MPI 4.0  

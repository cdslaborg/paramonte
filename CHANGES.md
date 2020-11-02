# ParaMonte kernel C/C++/Fortran release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases).  

## **Version 1.x.x**  

### **Version 1.5.0** -- Work in progress

**Enhancements**  

+   The ParaMonte C/C++/Fortran example build scripts `build.sh` and `build.bat` 
    now check for the existence of both of Intel and GNU compilers in the appropriate 
    order that is automatically inferred at the compilation time. Also, the dependencies 
    on the MPI compiler wrappers is now removed as the MPI libraries are not required to 
    build the ParaMonte examples, even in cases of parallel ParaMonte example builds.

### **Version 1.4.0** -- October 29, 2020

**Enhancements**  

+   The IO debugging info of all ParaMonte samplers have been enhanced. 
    In cases of wrong syntax or syntax-breaking input values in the simulation 
    output files, the error messages are now more informative and point directly 
    to the exact location of of error in the input file.  

+   The Integrated Autocorrelation (IAC) for sample refinement in ParaDRAM 
    sampler of ParaMonte is now set to the average of all variables' IAC values 
    instead of the maximum IAC value. This will lead to less aggressive decorrelation 
    of the final sample, which means significantly larger final sample sizes, without 
    compromising the i.i.d. property of the final refined sample. This behavior can 
    be reversed back to the original by specifying "max" or "maximum" along with 
    the requested refinement method, `SampleRefinementMethod = "batchmeans max"` 
    or `SampleRefinementMethod = "BatchMeans-max"` (case-insensitive).

### **Version 1.3.0** -- October 3, 2020

**Enhancements**  

+   A new simulation specification `overwriteRequested` has 
    been added to all ParaMonte samplers. If TRUE and the 
    ParaMonte sampler detects an existing set of old simulation 
    output files in the output path of the current simulation with 
    the same names as the output file names of the current simulation, 
    then, the ParaMonte sampler will overwrite the existing simulation files.  

**Compiler support**  

+   Intel Parallel Studio (>2018.0.0)  
+   GNU Compiler Collection (>7.0.0)  

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
        This can help make wiser decisions regarding the the number 
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

+   Intel Parallel Studio (>2018.0.0)  
+   GNU Compiler Collection (>7.0.0)  

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
+   The fully-deterministic restart functionality is now functional also when chainFileFormat="verbose" in ParaDRAM simulations.

**Compiler support**  

+   Intel Parallel Studio (>2018.0.0)  
+   GNU Compiler Collection (>7.0.0)  

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

+   Intel Parallel Studio (>2018.0.0)  
+   GNU Compiler Collection (>7.0.0)  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 8.3 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.166 20191121  
+   **macOS**: GNU 9.3 / Open-MPI 4.0  

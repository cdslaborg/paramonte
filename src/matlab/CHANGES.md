# ParaMonte MATLAB release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases) 
or [the ParaMonte page on MathWorks FileExchange central package repository](https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte).  

## **Version 3.x.x**  

### Version  3.0.0 -- Nov 13, 2024 (pre-release)  

**Major enhancements**  

+   This is a major release of the ParaMonte MATLAB library. 
    The library's usage syntax has changed in all programming language environments, including MATLAB. 

+   This update presents major performance, accuracy, and verification 
    enhancements to the ParaMonte kernel routines, in particular, 
    to the ParaDRAM sampler of the ParaMonte library.  

+   The ParaMonte MATLAB library now has a package structure, 
    allowing the use of different library functionalities without name clashes.

+   The ParaMonte MATLAB library now has a dedicated 
    [user and developer API documentation website](https://www.cdslab.org/paramonte/matlab/3/), 
    in addition to the generic documentation website containing the installation 
    and guidelines that apply to all programming language environments.

+   Numerous examples have been added to the library, all collected in 
    the `example` subfolder in the new ParaMonte MATLAB binary releases.

+   The visualization tools of the library have dramatically improved, expanded, and publicized.

+   The ParaMonte MATLAB samplers functionalities have significantly improved.

+   The ParaMonte MATLAB samplers can now take advantage of 
    MATLAB parallelism toolbox for shared-memory parallelism.
    For usage, see the sampling examples of the library.

+   The ParaMonte MATLAB MPI-parallel samplers are now significantly easier to use.
    For usage, see the sampling examples of the library.

+   Users are encouraged to test the most recent releases of the ParaMonte MATLAB available for 
    download at the repository's [GitHub release page](https://github.com/cdslaborg/paramonte/releases).

**Essential Dependency Compatibility**  

| Dependency                        | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| MATLAB >= R2023a (tested)         | ✅              | ✅            | ✅             | ✅            |  
| MATLAB <= R2022b (untested)       | ❓              | ❓            | ❓             | ❓            |  

**Optional Dependency Compatibility**  

| Dependency                        | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| Intel MPI (IMPI) >= 2021.8        | ✅              | ✅            | ✅             | ❌            |  
| MPICH MPI (MMPI) >= 3             | ❌              | ✅            | ✅             | ✅            |  
| OpenMPI (OMPI) >= 4               | ❌              | ✅            | ✅             | ✅            |  

## **Version 2.x.x**  

### Version  2.3.0 -- December 17, 2020  

**Major enhancements**  

+   This update presents several major performance, accuracy, 
    and verification enhancements to the ParaMonte kernel routines, 
    in particular, to the ParaDRAM sampler.  

+   An extensive set of over 866 tests have been added 
    that test all aspects of the ParaMonte kernel library.  

+   The issue of Windows file locking, that led to the occasional crashes of the 
    ParaDRAM and ParaDISE simulations in `multiChain` parallelism mode, is now resolved.  

+   The `ParaDRAM` class in `paramonte` is now also available 
    as `Paradram` and `paradram`, although the original label 
    will remain the default preferred method of ParaDRAM 
    object instantiation.  

### Version  2.2.1 -- November 15, 2020  

**Minor enhancements**  

+   Minor enhancements to the Kernel library 
    build scripts and dependencies management.  

+   More informative error messages are now printed 
    on MATLAB console if any error happens during the 
    ParaMonte MATLAB library setup on macOS for the first time.  

### Version  2.2.0 -- October 29, 2020  

**Enhancements**  

+   The `cmake` software dependency installation failure now 
    does not nullify the installation of other dependencies.

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

### Version  2.1.3 -- October 15, 2020  

**Minor enhancements**  

+   Further minor enhancements to the behavior of the 
    `checkForUpdate()` method of the `paramonte` class.  

### Version  2.1.2 -- October 15, 2020  

**Minor enhancements**  

+   The `checkForUpdate()` method of the `paramonte` 
    class now functions as expected.  

### Version  2.1.1 -- October 9, 2020  

**Minor enhancements**  

+   A Linux bug in the installation of the MPI library is now fixed.

### Version  2.1.0 -- October 3, 2020  

**Minor enhancements**  

+   A new simulation specification `overwriteRequested` has 
    been added to all ParaMonte samplers. If `True` and the 
    ParaMonte sampler detects an existing set of old simulation 
    output files in the output path of the current simulation with 
    the same names as the output file names of the current simulation, 
    then, the ParaMonte sampler will overwrite the existing simulation files.  

### Version  2.0.1 -- September 26, 2020  

**Minor enhancements**  

+   The guidelines for the installation of the 
    MPI library on macOS have been improved.  

+   The minor bug in GridPlot class method `rotateAxesLabels()` that caused 
    the `readSample()` , `readChain()`, `readChainMarkov()` to crash upon adding 
    Grid plots is now fixed.  

+   The minor bug in the naming of the ParaMonte kernel library files on macOS 
    (Darwin) is now fixed.  

### Version  2.0.0 -- September 22, 2020  

**Major enhancements to the ParaMonte / ParaDRAM sampler interfaces**  

+   The entire ParaMonte MATLAB interface library has been revamped.
    The new naming conventions, visualization, and computing tools 
    are significantly nicer to deal with and in some cases, orders 
    of magnitude faster than the previous major release.

+   The simulation output files reading is now completely overhauled. In particular, 
    the output file reader methods are now capable of handling input file paths that 
    point to a directory. In such cases, it will search the input directory for files 
    matching the requested file name pattern. If no input file is provided to the file 
    reader methods, the current working directory will be search for the potential 
    simulation files that match the requested pattern. 

+   Several new post-processing functionalities have now been added, such as
    the ability to seamlessly parse the contents of the output `*_report.txt`, 
    `*_restart.txt`, and `*_progress.txt` simulation files, in addition to the
    other output files (`*_sample.txt` and `*_chain.txt`) that could be parsed
    in the previous versions.

+   The newly-added `readRestart()` method is now added to the ParaDRAM sampler 
    class. User can now parse the contents of the output ASCII-format restart files. 
    This is particularly useful to visualize the dynamics of the ParaDRAM sampler class, 
    such as the evolution of the proposal distribution's location, shape, and covariance 
    matrix.  

+   The `GridPlot()` class now has two additional methods `setAxesLabels()` and 
    `setAxesLimits()` which can directly set the labels and limits of axes, hassle-free.

**Minor enhancements**  

+   The single value assignment to `spec.targetAcceptanceRate` component of a ParaDRAM object is now properly handles. 
    For example, the following code is valid as expected,  
    ```matlab  
    import paramonte as pm
    pmpd = pm.ParaDRAM()
    pmpd.spec.targetAcceptanceRate = 0.23 # this is now valid
    pmpd.spec.targetAcceptanceRate = [0.2, 0.3] # this is also valid, which limits the acceptance rate to the specified range
    ```  

+   The default background color in all plots is now `"white"`.  
+   The `rotateAxisLabels()` of the `GridPlot()` class is now renamed to `rotateAxesLabels()`.  

**Bug fixes**  

+   ParaDRAM `readChainMarkov()` no-output-option bug is now fixed. 
    When calling `readChainMarkov()`, user can now either provide the output variable or not.  

## **Version 1.x.x**  

### Version  1.1.0 -- June 5, 2020  

+   Enhancements and bug fixes to the kernel routines.  
+   Several major enhancements and bug fixes to the MATLAB kernel and interface routines.  
+   MatDRAM now supports fully-deterministic restart functionality.  

### Version  1.0.0 -- June 1, 2020 -- Initial release  

+   This is the first public release of the ParaMonte MATLAB library.  

**New features**  

+   ParaDRAM sampler: **Para**llel **D**elayed-**R**ejection **A**daptive Metropolis-Hastings **M**arkov Chain Monte Carlo Sampler.  
+   ParaMonte Interface to the MATLAB Programming language.  
+   ParaMonte simulation-output visualization via the ParaMonte MATLAB interface.  

# ParaMonte Python release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases) or [the ParaMonte page on the Python Package Index](https://pypi.org/project/paramonte/).  

## **Version 2.x.x**  

### Version  2.5.0 -- December 30, 2020  

**Major enhancements**  

+   This release is a a major step toward further portability 
    of the kernel routines of the `ParaMonte::Python` library. The kernel 
    library dependencies are now properly handled and recognized at runtime
    without such aggressive actions as permanently redefining the environmental
    path variables, most importantly, `PATH` and `LD_LIBRARY_PATH` on Linux/macOS.

+   The `ParaMonte::Python` library is now capable of recognizing the existing MPI
    libraries such as MPICH and OpenMPI on user's system and avoid further 
    installation of a new MPI library if it is deemed unnecessary.  

+   The ParaMonte kernel routines are now capable of handling user-input 
    file paths that contain white-space (blank) or other exotic characters.  

+   Enhancements to the `build()` function of the paramonte module.  

**Minor enhancements**  

+   Typo-fixes in the documentation of the library.

### Version  2.4.0 -- December 17, 2020  

**Major enhancements**  

+   This update presents several major performance, accuracy, 
    and verification enhancements to the ParaMonte kernel routines,
    in particular, to the ParaDRAM sampler.

+   An extensive set of over 866 tests have been added 
    that test all aspects of the ParaMonte kernel library.

+   The issue of Windows file locking, that led to the occasional crashes of the 
    ParaDRAM and ParaDISE simulations in `multiChain` parallelism mode, is now resolved.

+   All overflow / underflow exceptions are now properly handled.  

+   The `ParaDRAM` class in `paramonte` is now also available 
    as `Paradram` and `paradram`, although the original label 
    will remain the default preferred method of ParaDRAM 
    object instantiation.  

**Minor enhancements**  

+   The interactive mode for displaying plots is now automatically on. 
    The plots are not automatically displayed in `ipython` sessions.  

### Version  2.3.1 -- November 15, 2020  

**Minor enhancements**  

+   Minor enhancements to the Kernel library 
    build scripts and dependencies management.

### Version  2.3.0 -- October 29, 2020  

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
    the requested refinement method, `SampleRefinementMethod = "batchmeans-max"` 
    or `SampleRefinementMethod = "BatchMeans-maximum"` (case-insensitive).

### Version  2.2.5 -- October 20, 2020  

**Minor enhancements**  

+   Enhancements to the messages for incompatible architecture.  

### Version  2.2.4 -- October 18, 2020  

**Minor enhancements**  

+   The bug preventing the setting of 
    SpecDRAM specifications is now fixed.  

### Version  2.2.3 -- October 15, 2020  

**Minor enhancements**  

+   Further enhancements and corrections to the 
    `checkForUpdate()` function of `paramonte` module.  

### Version  2.2.2 -- October 14, 2020  

**Minor enhancements**  

+   The `checkForUpdate()` function of `paramonte` module now 
    works fine when a newer version of the software is available.  

### Version  2.2.1 -- October 11, 2020  

**Minor enhancements**  

+   Documentation enhancements and typo fixes.

### Version  2.2.0 -- October 11, 2020  

**Major enhancements**  

+   The `ParaMonte::Python` output file parsers method are now 
    capable of parsing simulation output contents directly from 
    the web. All that is needed, is to provide the link to the 
    web file as the input file name to the simulation output 
    file parser methods (e.g., `readSample()`, `readChain()`, 
    `readReport()`, `readProgress()`, `readRestart()`, ...).

**Minor enhancements**  

+   A bug in the y-axis variable names in the heatmap plot is now fixed.

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

### Version  2.0.9 -- October 2, 2020  

**Minor enhancements**  

+   Minor correction to the value of `__version__`, 
    now representing solely the version number.  

+   A simple example-usage Python script is now 
    added to the README.md file of the package.  

### Version  2.0.8 -- September 29, 2020  

**Minor enhancements**  

+   Enhanced error messages for situations when 
    the MPI library cannot be found on the system.  

### Version  2.0.7 -- September 26, 2020  

**Minor enhancements**  

+   The guidelines for the installation of the MPI 
    library on macOS and Linux have been improved.  

### Version  2.0.6 -- September 25, 2020  

**Minor enhancements**  

+   The explicit dependencies on `scipy`, `matplotlib`, and `seaborn` are 
    now removed from the PyPI setup file of the ParaMonte library as these 
    are only required for the post-processing and visualizations of the 
    simulation results. From now on, only `numpy` and `pandas` are the 
    minimally-required Python modules, and practically, only `numpy`.  

+   Two new functions `verifyDependencyVersion()` and `getDependencyVersion()`
    are now added to the library that can check for the existence of the 
    ParaMonte library's visualization dependencies and their required 
    minimum versions.  

+   The `seaborn` Python library has now decided to deprecate the `distplot()` 
    function. The corresponding visualization method in the ParaMonte library 
    has been now updated to a more appropriate name and underlying function.  

### Version  2.0.4 -- September 22, 2020  

**Minor enhancements**  

+   The output of the plotting functions is now stored as a list in 
    the `currentFig` temporary component of the visualization objects.
    This way, access to multiple individual objects on the active plot 
    is maintained instead of only the last object. Overall, this is a 
    minor change that will not cause any noticeable change in the 
    behavior of the library in almost in all use cases.

+   A minor bug regarding the input value for the `outputDelimiter` 
    attribute of the `spec` component of the `ParaMonteSampler()` class,  
    used in the `readTabular()` internal method, is now fixed.

### Version  2.0.3 -- September 11, 2020  

**Minor enhancements**  

+   Minor enhancement to `checkForUpdate()` method of 
    the `paramonte` module.  

### Version  2.0.2 -- September 11, 2020  

**Minor enhancements**  

+   Minor enhancement to `checkForUpdate()` method of 
    the `paramonte` module.  

### Version  2.0.1 -- September 10, 2020  

**Minor enhancements**  

+   LGPL3 LICENSE is now switched to MIT LICENSE.md file.  

+   A fix to the `brew` software installation now avoids 
    the seemingly-unavoidable crash.  

### Version  2.0.0 -- September 6, 2020  

**Major enhancements to the ParaMonte / ParaDRAM sampler interfaces**  

+   The entire ParaMonte Python interface library has been revamped.
    The new naming conventions, visualization, and computing tools 
    are significantly nicer to deal with and in some cases, orders 
    of magnitude faster than the previous major release.

+   The kernel density estimates and visualization tools are now on average 
    **100 times or more faster than the previous release of the library**.

+   Several new post-processing functionalities have now been added, such as
    the ability to seamlessly parse the contents of the output `*_report.txt`, 
    `*_restart.txt`, and `*_progress.txt` simulation files, in addition to the
    other output files (`*_sample.txt` and `*_chain.txt`) that could be parsed
    in the previous versions.

+   The new major release also includes 3D visualization tools, such as 3D 
    line, scatter, or line+scatter plots as well as fast 2D and 3D kernel 
    density estimate contour plotting tools.

**Minor enhancements**  

+   The simulation output files reading is now completely overhauled. In particular, 
    the output file reader methods are now capable of handling input file paths that 
    point to a directory. In such cases, it will search the input directory for files 
    matching the requested file name pattern. If no input file is provided to the file 
    reader methods, the current working directory will be search for the the potential 
    simulation files that match the requested pattern. 

+   The error-signaling behavior of the library now is very much controlled, that is, 
    upon code failure, it does not automatically shutdown the Python kernel in Jupyter 
    Notebooks. The library now simply throws an error message upon failing instead of 
    restarting the environment.  

+   The single value assignment to `spec.targetAcceptanceRate` component of a ParaDRAM 
    object is now properly handled. For example, the following code is valid as expected,  
    ```python  
    import paramonte as pm
    pmpd = pm.ParaDRAM()
    pmpd.spec.targetAcceptanceRate = 0.23 # this is now valid
    pmpd.spec.targetAcceptanceRate = [0.2, 0.3] # this is also valid, which limits the acceptance rate to the specified range
    ```  

+   The minimum required dependency versions are now raised to the following,  
    ```python  
    python_requires = ">=3.5"
    install_requires = [ "numpy>=1.18.0"
                       , "scipy>=1.4.0"
                       , "pandas>=1.0.0"
                       , "seaborn>=0.10.0"
                       , "matplotlib>=3.2.0"
                       ]
    ```  


## **Version 1.x.x**  

### Version  1.1.1 -- June 7, 2020  

**Minor enhancements**  

+   The `_ScatterLinePlot` dangling class is removed from the package.  

### Version  1.1.0 -- June 1, 2020  

+   Major enhancements to the ParaMonte kernel library.  
+   Major bug fixes in the ParaMonte Python library.  
+   The ParaMonte kernel and Python interface versions are now reposted separately as components of the paramonte module. 

### Version  1.0.12 -- April 6, 2020  

+   Minor enhancements and bug fixes to the kernel routines.

### Version  1.0.11 -- April 4, 2020  

+   Minor enhancements and bug fixes to the GridPlot.

### Version  1.0.10 -- March 28, 2020  

+   Minor bug fix.

### Version  1.0.9 -- March 27, 2020  

+   Minor enhancements.

### Version  1.0.8 -- March 27, 2020  

+   Minor enhancements.

### Version  1.0.7 -- March 26, 2020  

+   Minor corrections.

### Version  1.0.6 -- March 22, 2020  

+   Minor bug fixes.

### Version  1.0.5 -- March 21, 2020  

+   Minor bug fix.

### Version  1.0.4 -- March 20, 2020  

+   support for macOS (Darwin) added.

### Version  1.0.3 -- February 13, 2020  

+   Minor bug fixes.

### Version  1.0.2 -- February 13, 2020  

+   Minor bug fixes to the parallel routines.

### Version  1.0.1 -- February 13, 2020  

+   Minor bug fixes.

### Version  1.0.0 -- January 1, 2020 -- Initial release  

+   This is the first public release of the ParaMonte library.  

**New features**  
+   ParaDRAM sampler: **Para**llel **D**elayed-**R**ejection **A**daptive Metropolis-Hastings **M**arkov Chain Monte Carlo Sampler.  
+   ParaMonte Interface to the Python Programming languages.  
+   ParaMonte simulation-output visualization via the ParaMonte Python interface.  

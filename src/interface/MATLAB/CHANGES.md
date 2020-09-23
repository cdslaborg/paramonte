# ParaMonte MATLAB release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases).  

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
    reader methods, the current working directory will be search for the the potential 
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

+   ParaDRAM `readMarkovChain()` no-output-option bug is now fixed. 
    When calling `readMarkovChain()`, user can now either provide the output variable or not.  

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

# ParaMonte Python release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases) or [the ParaMonte page on the Python Package Index](https://pypi.org/project/paramonte/).  

## **Version 1.x.x**  

### Version  1.1.2 -- Work in Progress

**Minor enhancements**  

+   The error-signaling behavior of the library now is very much controlled, that is, 
    upon code failure, it does not automatically shutdown the Python kernel in Jupyter Notebooks. 
    The library now simply throws an error message upon failing instead of restarting the environment.  

+   The single value assignment to `spec.targetAcceptanceRate` component of a ParaDRAM object is now properly handles. 
    For example, the following code is valid as expected,  
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

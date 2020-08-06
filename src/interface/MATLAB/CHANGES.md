# ParaMonte MATLAB release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases).  

## **Version 1.x.x**  

### Version  1.1.1 -- Work in Progress  

**bug fixes**  
+   ParaDRAM `readMarkovChain()` no-output-option bug is now fixed. User can now either provide the output variable or not when calling `readMarkovChain()`.  

**New features**  
+   new `readReport()` method now added to the ParaDRAM sampler class. User can now parse the contents of the output report file.  

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

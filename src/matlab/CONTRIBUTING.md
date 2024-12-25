
Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started with **contributing to the ParaMonte MATLAB library**.  

## Initial Steps [⛓](#initial-steps-)

+   First, read the [general development guidelines](../../CONTRIBUTING.md). 
+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. 
    Make sure you find an open issue **about the MATLAB routines** and that you do not duplicate someone else's work.  
+   If your contribution does not exist as an issue, post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement/feature request(s) or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  

## ParaMonte MATLAB Development Conventions [⛓⛓](#paramonte-matlab-development-conventions-)

+   Never use any MATLAB feature that is younger than 3 years.  
    **Why?** The ParaMonte MATLAB library strives to be as compatible 
    as possible with all MATLAB releases younger than five years.  

+   Always follow the Doxygen-for-MATLAB documentation rules for documenting the library routines.  
    You can always copy and paste a template from an existing procedure documentation to create new.  
    **Why?** The ParaMonte MATLAB library uses an extension of [Doxygen for MATLAB](https://github.com/simgunz/doxymatlab)
    for documentation purposes.

+   Always use MATLAB packages (corresponding to `+folder` structures) to organize and structure the newly added functionalities within the library.
    Anything that must not be exposed to the end users, must not appear in a package folder directly. 
    Instead, it must be placed in a regular folder to remain hidden from the end users.
    **Why?** MATLAB lacks the concept of modules. Instead, it relies on the concept of
    packages organized as a nested structure of folders whose names begin with a `+`.

## Final Steps [⛓⛓](#final-steps-)

Once you have implemented your contributions,  

+   Do not forget to test your contributions by adding new tests to the library's unit-testing framework.  
+   Also, generate a code coverage report to ensure your contributions do not lower the overall code coverage of the library routines.  
+   Follow the [generic contribution guidelines](../../CONTRIBUTING.md/#all-contributors) to submit and merge your contributions with the library's main branch on GitHub. 

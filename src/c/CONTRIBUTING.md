
Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started with **contributing to the ParaMonte C library**.  

## Initial Steps [⛓](#initial-steps-)

+   First, read the [general development guidelines](../../CONTRIBUTING.md). 
+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. 
    Make sure you find an open issue **about the C routines** and that you do not duplicate someone else's work.  
+   If your contribution does not exist as an issue, post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement/feature request(s) or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  

## ParaMonte C Development Conventions [⛓⛓](#paramonte-c-development-conventions-)

Pay careful attention to the following conventions used in developing the C routines.

### Preprocessor Directives [⛓⛓⛓](#preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#preprocessor-directives-).

#### Platform Preprocessor Directives [⛓⛓⛓⛓](#platform-preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#platform-preprocessor-directives-).

#### Compiler Vendor Preprocessor Directives [⛓⛓⛓⛓](#compiler-vendor-preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#compiler-vendor-preprocessor-directives-).

#### Parallelism Preprocessor Directives [⛓⛓⛓⛓](#parallelism-preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#parallelism-preprocessor-directives-).

#### Library-Type Preprocessor Directives [⛓⛓⛓⛓](#library-type-preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#library-type-preprocessor-directives-).

#### Library-Build Preprocessor Directives [⛓⛓⛓⛓](#library-build-preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#library-build-preprocessor-directives-).

#### Library-Interface Preprocessor Directives [⛓⛓⛓⛓](#library-interface-preprocessor-directives-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#library-interface-preprocessor-directives-).

### Coding Style Conventions [⛓⛓⛓](#coding-style-conventions-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#coding-style-conventions-).

### Controlling the Runtime Checks [⛓⛓⛓](#controlling-the-runtime-checks-)

Follow the rules specified for the [ParaMonte Fortran development](../fortran/CONTRIBUTING.md#controlling-the-runtime-checks-).

## Final Steps [⛓⛓](#final-steps-)

Once you have implemented your contributions,  

+   Do not forget to test your contributions by adding new tests to the library's unit-testing framework.  
+   Also, generate a code coverage report to ensure your contributions do not lower the overall code coverage of the library routines.  
+   Follow the [generic contribution guidelines](../../CONTRIBUTING.md/#all-contributors) to submit and merge your contributions with the library's main branch on GitHub. 

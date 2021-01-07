
Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started with **contributing to the ParaMonte kernel routines**.  

## Initial steps  

+   First, read the [general development guidelines](CONTRIBUTING.md). 
+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. 
    Make sure you find an open issue **about the kernel routines** and that you do not duplicate someone else's work.  
+   If your contribution does not exist as an issue, post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement/feature request(s), or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  

## Kernel development conventions  

Pay careful attention to the following conventions used in the development of the kernel routines.

### Preprocessor directives  

The kernel routines of ParaMonte heavily rely on the compiler preprocessing directives to implement multiple completely different parallelism paradigms (serial, MPI, Coarray) that are specifically tailored for the needs of multiple completely different programming languages (C, C++, Fortran, Julia, MATLAB, Python, R, ...) for different builds (`debug`, `testing`, `release`), for different operating systems (`Windows`, `Linux`, `macOS`, `WSL`), and many more.  

Despite significant efforts to minimize the use of preprocessor directive in the kernel routines, their frequent usage is unavoidable as they greatly minimize the cost of development. Simultaneously, however, the preprocessor directive can increase confusion and and lead to implicit bugs that are hard to detect.  

> **WARNING**  
> Any code section that requires some specific operating system, platform, compiler, programming language, compiler settings, or parallelism paradigms must be fenced with the appropriate preprocessor macros such that it is executed only in the appropriate circumstances. Failure to do so, will lead to a non-portable codebase. In general, there must exist an `#else` preprocessor section for any fenced code section that appropriately handles other possible scenarios, if the other scenarios are not possible, it is followed by an `#error \"error message\"` preprocessor directive.  

> **IMPORTANT**  
> The preprocessor macros must be always defined with uppercase characters. This is to easily distinguish macros from other variables in the regular sections of the code.  

> **IMPORTANT**  
> Define macros as propositions that evaluate to either TRUE or FALSE. For example, `OS_IS_WIDOWS` or `MPI_ENABLED`.  

> **TIP**  
> When developing or debugging a particular feature that is bound to a specific set of preprocessor macros, you can request the compiler to generate intermediate preprocessed files that are free from preprocessor directives and let you focus on the problem better. In doing so, keep in mind that every change you make to the intermediate codes has be carefully implemented in the original set of codes.  

#### Platform preprocessor directives  

If there is any section of the code that must be exclusively executed on a particular platform or Operating System (OS), it must be properly fenced by the appropriate preprocessor flag.

1.  The `OS_IS_WINDOWS` preprocessor macro must be defined for any code section that exclusively belongs to Windows operating systems.  
1.  The `OS_IS_DARWIN` preprocessor macro must be defined for any code section that exclusively belongs to macOS operating systems.  
1.  The `OS_IS_LINUX` preprocessor macro must be defined for any code section that exclusively belongs to Linux operating systems.  
1.  The `OS_IS_WSL` preprocessor macro must be defined for any code section that exclusively belongs to Microsoft Subsystem for Linux.  
    > **NOTE**  
    > The `OS_IS_WSL` preprocessor macro is mostly useful for defining code sections that are written for compatibility with Microsoft Subsystem for Linux Version 1 (**WSL1**). For example, WSL1 is known to have issues with passing internal procedures to external routines since [WSL has a non-executable stack](https://github.com/Microsoft/WSL/issues/3083). Although, [this issue](https://github.com/Microsoft/WSL/issues/286) appears to have been resolved in [WSL2](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux#WSL_2), as a developer, you should be mindful of all users on all platforms.  

Several issues also need special care when developing ParaMonte on Windows OS:  

+   On Windows systems, **file-locking** by a processor often leads to problems and program crashes in parallel simulations. Since, currently only Intel compilers are supported and tested with the ParaMonte library, a quick remedy is activate the shared file manipulation via the Intel compiler `SHARED` argument that can be passed to `open()` statements. For example,  
    ```fortran  
    open( unit      = self%SampleFile%unit            &
        , file      = self%SampleFile%Path%original   &
        , status    = self%SampleFile%status          &
        , iostat    = self%SampleFile%Err%stat        &
    #if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
        , SHARED                                      &
    #endif
        , position  = self%SampleFile%Position%value  )
    ```  
+   To properly build shared (DLL) libraries on Windows OS, every procedure name must be correctly extracted via Intel compiler `DLLEXPORT` directive. Currently, all names are manually extracted immediately below each procedure's interface. This is normally not an issue, as long the symbols are correctly exported when the procedure is being developed. However, this can quickly become a challenge if the developer uses a wrong symbol name for the procedure or completely forgets to export the symbol.  
    > **NOTE**  
    > A potential solution to this problem is to enable builds with CMAKE on Windows and [request CMAKE to automatically export all symbols](https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/). Indeed, extending the current cmake build system of the ParaMonte library to Windows would be a great contribution to the package, if you can accomplish it. Before doing so, check the issues page of the project to make sure you would not duplicate someone else's work.  

#### Compiler suite preprocessor directives  

Generally, one should avoid the use of code that bind the library to a particular compiler suite. Sometimes, however, this is unavoidable. Currently, the library supports builds with the **Intel** and **GNU** compiler suites.

1.  The `INTEL_COMPILER_ENABLED` preprocessor macro must be defined for any code section that requires the Intel compiler to recognize the syntax. This preprocessor flag is automatically defined by the ParaMonte build scripts and passed to the Intel preprocessor when the Intel compiler is used to build the library.  
1.  The `GNU_COMPILER_ENABLED` preprocessor macro must be defined for any code section that requires the GNU compiler to recognize the syntax. This preprocessor flag is automatically defined by the GNU preprocessor on all platforms.  

#### Parallelism preprocessor directives  

There are currently two preprocessor directives that determine the type of parallelism to which the code section sandwiched by the preprocessor macro belongs.  

1.  The `MPI_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to MPI-parallelized versions of the kernel routines.  
1.  The `CAF_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to CAF-parallelized versions of the kernel routines (Coarray parallelism).  
1.  The `OMP_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to OMP-parallelized versions of the kernel routines (OpenMP parallelism). As of Dec 2020, there is no section of the kernel routines uses OpenMP parallelism.  

For example, if a section is MPI-parallelized, it must be fenced with the following preprocessor directive,  

```fortran  
#if defined MPI_ENABLED
call mpi_initialized( isInitialized, ierrMPI )
#endif
```  

or, if a section must be executed in either MPI or Coarray parallelisms, it must be fenced via the following preprocessor directive,  

```fortran  
#if defined CAF_ENABLED || defined MPI_ENABLED
integer :: imageID
#endif
```  

#### Library type preprocessor directives  

Any code section that must be activated for a particular library build type (`static` vs. `shared`) must be fenced with appropriate preprocessing macros.  

1.  The `DLL_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to shared library builds. Although `DLL` is reminiscent of Windows shared library files, the ParaMonte build scripts define this preprocessor macro for all shared library builds irrespective of the platform (Windows, Linux, macOS, ...).  
    > **IMPORTANT**  
    > Use this preprocessor macro in combination with **`IFORT_ENABLED`** to enable DLL symbol export by the Intel **`ifort`** compiler on macOS and Windows operating systems.  

#### Library build preprocessor directives  

Occasionally, it is necessary to define sections of code that should run only in a particular library build. Currently, the ParaMonte library build scripts support three different library builds (`debug`, `testing`, and `release`). Any code section that does not belong to the release (production) build must be appropriately fenced with appropriate preprocessor flags.  

1.  The `DEBUG_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the `debug` library builds. Use this macro to fence sections of code that help with the debugging of the code.
1.  The `TESTING_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the `testing` library builds. Use this macro to fence sections of code that must be activated during the testing of the library.
    > **TIP**  
    > The `TESTING_ENABLED` preprocessor macro often appears together with `DEBUG_ENABLED` and `CODECOV_ENABLED`. The testing mode deactivates all compiler optimizations, but also prevents the addition of the debugging symbols to the object files. This leads to faster compilation of the source files which is often desired for testing purposes.  

#### Library interface preprocessor directives  

To make the ParaMonte an inter-operable cross-language library, it is necessary to tailor sections of code for the particular languages of interest. The kernel connection to all languages other than Fortran is provided via the `FCI_ENABLED` preprocessor macro in the ParaMonte library.  

1.  The `CFI_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible from any programming language other than Fortran, most importantly, the C programming language.  

In addition, some code sections might have to be exclusively executed when the library is invoked from a particular programming language. In such cases, the relevant code section must be properly fenced with the appropriate language-specific preprocessor macros.  

1.  The `C_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the C programming language.  
1.  The `CPP_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the C++ programming language.  
1.  The `FORTRAN_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Fortran programming language.  
1.  The `JULIA_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Julia programming language.  
1.  The `MATLAB_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the MATLAB programming language.  
1.  The `MATHEMATICA_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Mathematica programming language.  
1.  The `PYTHON_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Python programming language.  
1.  The `R_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the R programming language.  

### Coding style conventions  

The following coding style are enforced within the kernel files. 
If you do not follows these rules in your contribution, please provide an explanation and justification of the alternative approach that you have taken in your contribution.  

+   [Strict naming conventions are enforced](CONTRIBUTING.md/#general-coding-style-conventions) within the entire library, including the kernel routines.

+   Strict semantic compliance with the latest Fortran standard (2018).

+   Strict source-code compliance with the latest Fortran standard. However, there are important vital exceptions to this rule,  
    +   It often preferable and sometimes essential to continue the source code lines beyond the current maximum line length limit specified by the standard, which is 132 characters. This maximum line length limit restriction is planned to be lifted in the next release of Fortran 202X standard. Nevertheless, a line length of approximately 200 characters is about where you should seriously think of breaking and continuing the code on multiple lines, if needed.  
    +   Strict parallelism compliance with the Coarray Fortran parallelism is not required (via coarrays). In other words, you can implement parallelism that is provided by external libraries, the most prominent example of which is the Message Passing Interface (MPI) Library. However, in doing so, you will have properly fence the MPI section with the proper preprocessor macro as explained in the previous section.  
    +   The minimal and efficient use of preprocessor macros and directives is allowed and encouraged to avoid code duplication and improve productivity.  

+   All Coarray variables or objects (that are communicated between processes) must begin with the `co_` prefix. Example: `co_logFunc`.  
    **Why?** This helps easily identify object that are responsible for inter-process communications in the Coarray implementation of the library.  

+   All **m**odule **v**ariables or module objects must begin with `mv_`, unless there is a good justification for not doing so. Example: `mv_state`, `comv_logFuncState`.  
    Why? module variables are global, at least within the scope of the module, and at most, everywhere in the program. Global variables are extremely susceptible of creating extremely difficult hard-to-identify bugs in the code and are hard to follow. Adding the prefix `mv_` to all global variables, resolves at least part of the problem by explicitly declaring the variable as global to the developer via its name.  

+   All **m**odule **c**onstant variables that are supposed to be defined only once throughout the entire runtime must be prefixed with `mc_` instead of `mv_`. Example: `mc_imageCount = 3`.  
    **Why?** Adding the prefix `mc_` helps other developers to realize that the object is not supposed to change value at runtime beyond the first initialization. Note that this type of variable is different from a compile-time constant, since its value must be defined at runtime, which could be different from one run to another. However, once set, the value is expected to remain constant throughout the current simulation.  

+   All constants (parameters) are upper-case separated by underscore. Example: `FILE_EXT = ".txt"`.  

+   All module names must begin with an uppercase letter and more importantly, must end with the prefix `_mod`. Example: `ParaMCMC_mod`.  
    **Why?** This helps create multiple different entity names from a single base name. For example,  
    ```fortran  
    use String_mod, only: String_type
    type(String_type) :: string
    ```  

+   Each module must be stored in a separate source file and the name of the source file must be the same as the module name, unless there is a good justification to do it otherwise. This is a strict rule to maintain the integrity of the library source file and the build scripts.  

+   All submodule names must begin with an uppercase letter and more importantly, must end with the prefix `_smod`. Example: `Routines_smod`.  

+   If a module has submodules, the name of the submodule file should follow the module name after `@`. For example the following submodule,  
    ```fortran  
    submodule (Decoration_mod) Routines_smod
    ...
    end submodule Routines_smod
    ```  
    will have to be stored in a source file with the name `Decoration_mod@Routines_smod.f90`.  

+   All type and class names must begin with an uppercase letter and more importantly, must end with the prefix `_type`. Example: `ParaMCMC_type`.  
    **Why?** This helps create multiple different entity names from a single base name. For example,  
    ```fortran  
    use String_mod, only: String_type
    type(String_type) :: string
    ```  

## Final steps  

Once you have implemented your contributions,  

+   Do not forget to test your contributions by adding new kernel tests to the unit-testing framework of the library.  
+   Also, generate code coverage report to ensure your contributions do not lower the overall code coverage of the kernel routines.  
+   Follow the [generic contribution guidelines](CONTRIBUTING.md/#all-contributors) to submit and merge your contributions with the main branch of the library on GitHub. 

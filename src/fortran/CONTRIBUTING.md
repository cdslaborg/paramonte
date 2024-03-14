
Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started with **contributing to the ParaMonte Fortran library**.  

## Initial Steps [⛓](#initial-steps-)

+   First, read the [general development guidelines](../../CONTRIBUTING.md). 
+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. 
    Make sure you find an open issue **about the Fortran routines** and that you do not duplicate someone else's work.  
+   If your contribution does not exist as an issue, post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement/feature request(s) or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  

## ParaMonte Fortran Development Conventions [⛓⛓](#paramonte-fortran-development-conventions-)

Pay careful attention to the following conventions used in developing the Fortran routines.

### Preprocessor Directives [⛓⛓⛓](#preprocessor-directives-)

The Fortran routines of ParaMonte heavily rely on the compiler preprocessing directives to implement multiple completely different parallelism paradigms (serial, MPI, Coarray) that are specifically tailored for the needs of multiple completely different programming languages (C, C++, Fortran, Julia, MATLAB, Python, R, ...) for different builds (`debug`, `testing`, `release`), for different operating systems (`Windows`, `Linux`, `macOS`, `WSL`), and many more.  

Despite significant efforts to minimize preprocessor directives in the Fortran routines, their frequent usage is unavoidable as they greatly minimize the cost of development. 
Simultaneously, however, the preprocessor directive can increase confusion and lead to implicit bugs that are hard to detect.  

>   **WARNING**  
>   Any code section that requires some specific operating system, platform, compiler, programming language, compiler settings, or parallelism paradigms must be fenced with the appropriate preprocessor macros such that it is executed only in the appropriate circumstances. Failure to do so will lead to a non-portable codebase. In general, an `#else` preprocessor section must exist for any fenced code section that appropriately handles other possible scenarios if the other scenarios are not possible, followed by an `#error \"error message\"` preprocessor directive.  

>   **WARNING**  
>   Macro names, if they are to be used as conditions, must be always given a value of `1` if true.<br>
>   This is specially important if the macro is defined from within the source code.<br>
>   **Why?** Such value assignment eliminates reliance of the conditional construct on the *definition*
>   status of the macro as the condition. For example, the following conditional construct executes
>   because `CHECK_ENABLED` is assigned a value of `1`. Example:<br>
>   ```fortran  
>   #define CHECK_ENABLED 1  
>   #if CHECK_ENABLED  
>       if (.not. present(err) .and. choDia(1) < 0._RK) error stop  
>   #endif  
>   ```  

>   **IMPORTANT**  
>   The preprocessor macros must always be defined with uppercase characters. This is to easily distinguish macros from other variables in the regular sections of the code.  

>   **IMPORTANT**  
>   Define macros as propositions that evaluate to either TRUE or FALSE. For example, `WINDOWS_ENABLED` or `MPI_ENABLED`.  

>   **TIP**  
>   When developing or debugging a particular feature bound to a specific set of preprocessor macros, you can request the 
>   compiler to generate intermediate preprocessed files free from preprocessor directives and let you focus on the problem better. 
>   In doing so, remember that every change you make to the intermediate codes must be carefully implemented in the original set of codes.  

#### Platform Preprocessor Directives [⛓⛓⛓⛓](#platform-preprocessor-directives-)

Suppose a code section must be exclusively executed on a particular platform or Operating System (OS). In that case, it must be appropriately fenced by the appropriate preprocessor flag.

1.  The `WINDOWS_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the Windows operating system.  
1.  The `DARWIN_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to macOS operating systems.  
1.  The `LINUX_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to Linux operating systems.  
1.  The `WSL_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to Microsoft Subsystem for Linux.  
    > **NOTE**  
    > The `WSL_ENABLED` preprocessor macro is mostly useful for defining code sections written for compatibility with Microsoft Subsystem for Linux Version 1 (**WSL1**). For example, WSL1 is known to have issues with passing internal procedures to external routines since [WSL has a non-executable stack](https://github.com/Microsoft/WSL/issues/3083). Although [this issue](https://github.com/Microsoft/WSL/issues/286) appears to have been resolved in [WSL2](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux#WSL_2), as a developer, you should be mindful of all users on all platforms.  

Several issues also need special care when developing ParaMonte on Windows OS:  

+   On Windows systems, **file-locking** by a processor often leads to problems and program crashes in parallel simulations. Since currently, only Intel compilers are supported and tested with the ParaMonte library, a quick remedy is to activate the shared file manipulation via the Intel compiler `SHARED` argument that can be passed to `open()` statements. For example,  
    ```fortran  
    open( unit      = self%sampleFile%unit            &
        , file      = self%sampleFile%Path%original   &
        , status    = self%sampleFile%status          &
        , iostat    = self%sampleFile%err%stat        &
    #if INTEL_ENABLED && WINDOWS_ENABLED
        , SHARED                                      &
    #endif
        , position  = self%sampleFile%Position%val  )
    ```  
+   To build shared (DLL) libraries on Windows OS properly, every procedure name must be correctly extracted via the Intel compiler `DLLEXPORT` directive. Currently, all names are manually extracted immediately below each procedure's interface. This is usually not an issue as long the symbols are correctly exported when the procedure is being developed. However, this can quickly become a challenge if the developer uses wrong symbol name for the procedure or completely forgets to export the symbol.  
    > **NOTE**  
    > A potential solution to this problem is to enable builds with CMake on Windows and [request CMAKE to export all symbols automatically](https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/). Indeed, extending the current CMake build system of the ParaMonte library to Windows would be a great contribution to the package, if you can accomplish it. Before doing so, check the project's issues page to ensure you do not duplicate someone else's work.  

#### Compiler Vendor Preprocessor Directives [⛓⛓⛓⛓](#compiler-vendor-preprocessor-directives-)

Generally, one should avoid using code that binds the library to a particular compiler suite.<br>
Sometimes, however, this is unavoidable. Currently, the library supports builds with the **Intel** and **GNU** compiler suites.<br>

1.  The `INTEL_ENABLED` preprocessor macro must be defined for any code section that requires the Intel compiler to recognize the syntax. This preprocessor flag is automatically defined by the ParaMonte build scripts and passed to the Intel preprocessor when the Intel compiler is used to build the library.  
1.  The `GNU_ENABLED` preprocessor macro must be defined for any code section that requires the GNU compiler to recognize the syntax. This preprocessor flag is automatically defined by the GNU preprocessor on all platforms.  

#### Parallelism Preprocessor Directives [⛓⛓⛓⛓](#parallelism-preprocessor-directives-)

Currently, two preprocessor directives determine the type of parallelism to which the code section sandwiched by the preprocessor macro belongs.  

1.  The `MPI_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to MPI-parallelized versions of the Fortran routines.  
1.  The `CAF_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to CAF-parallelized versions of the Fortran routines (Coarray parallelism).  
1.  The `OMP_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to OMP-parallelized versions of the Fortran routines (OpenMP parallelism).    

For example, if a section is MPI-parallelized, it must be fenced with the following preprocessor directive,  

```fortran  
#if MPI_ENABLED
call mpi_initialized(isInitialized, ierrMPI)
#endif
```  

or, if a section must be executed in either MPI or Coarray parallelisms, it must be fenced via the following preprocessor directive,  

```fortran  
#if CAF_ENABLED || MPI_ENABLED
integer :: imageID
#endif
```  

#### Library-Type Preprocessor Directives [⛓⛓⛓⛓](#library-type-preprocessor-directives-)

Any code section that must be activated for a particular library build type (`static` vs. `shared`) must be fenced with appropriate preprocessing macros.  

1.  The `DLL_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to shared library builds. 
    Although `DLL` is reminiscent of Windows shared library files, the ParaMonte build scripts define this preprocessor macro for all shared library builds irrespective of the platform (Windows, Linux, macOS, ...).  
    > **IMPORTANT**  
    > Use this preprocessor macro in combination with **`INTEL_ENABLED`** to enable DLL symbol export by the Intel **`ifort`** compiler on macOS and Windows operating systems.  

#### Library-Build Preprocessor Directives [⛓⛓⛓⛓](#library-build-preprocessor-directives-)

Occasionally, it is necessary to define code sections that should be included only for a particular library build. 
Currently, the ParaMonte library build scripts support three different library builds (`debug`, `testing`, and `release`). 
Any code section not belonging to the release (production) build must be appropriately fenced with appropriate preprocessor flags.  

1.  The `DEBUG_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the `debug` library builds. 
    Use this macro to fence sections of code that help debug the code.
1.  The `TESTING_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the `testing` library builds. 
    Use this macro to fence sections of code that must be activated during the library testing.
    > **TIP**  
    > The `TESTING_ENABLED` preprocessor macro often appears together with `DEBUG_ENABLED` and `CODECOV_ENABLED`. 
    > The testing mode deactivates all compiler optimizations but also prevents the addition of the debugging symbols to the object files. 
    > This leads to faster compilation of the source files, which is often desired for testing purposes.  

#### Library-Interface Preprocessor Directives [⛓⛓⛓⛓](#library-interface-preprocessor-directives-)

Making ParaMonte an interoperable cross-language library requires customizing code sections for the particular languages of interest. 
The Fortran connection to all languages other than Fortran is provided via the `FCI_ENABLED` preprocessor macro in the ParaMonte library.  

1.  The `CFI_ENABLED` preprocessor macro must be defined for any code section accessible from any programming language other than Fortran, most importantly, the C programming language.  

In addition, some code sections must be exclusively executed when the library is invoked from a particular programming language. 
In such cases, the relevant code section must be properly fenced with the appropriate language-specific preprocessor macros.  

1.  The `C_ENABLED` preprocessor macro must be defined for any code section to be accessible only when the library is built for the C programming language.  
1.  The `CPP_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the C++ programming language.  
1.  The `FORTRAN_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the Fortran programming language.  
1.  The `JULIA_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the Julia programming language.  
1.  The `MATLAB_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the MATLAB programming language.  
1.  The `MATHEMATICA_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the Mathematica programming language.  
1.  The `PYTHON_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the Python programming language.  
1.  The `R_ENABLED` preprocessor macro must be defined for any code section accessible only when the library is built for the R programming language.  

### Coding Style Conventions [⛓⛓⛓](#coding-style-conventions-)

The following coding styles are enforced within the Fortran files. 
If you do not follow these rules in your contribution, please explain and justify the alternative approach you have taken in your contribution.  

+   Always prefer clarity to conciseness, specially when choosing variable, argument, or procedure names.

+   [Strict naming conventions are enforced](../../CONTRIBUTING.md/#general-coding-style-conventions) within the entire library, including the Fortran routines.

+   Strict semantic compliance with the latest Fortran standard (2018).

+   Strict source-code compliance with the latest Fortran standard. However, there are important exceptions to this rule,  
    +   It is often preferable and sometimes essential to continue the source code lines beyond the current maximum line length limit specified by the standard, which is 132 characters. 
        This maximum line length limit restriction is lifted in the most recent release of the Fortran standard (2023). 
        Nevertheless, a line length of approximately 200-300 characters is where you should seriously consider breaking and continuing the code on multiple lines if needed.  
    +   Strict parallelism compliance with the Coarray Fortran parallelism is not required (via Coarrays). 
        In other words, you can implement parallelism provided by external libraries, the most prominent example of which is the Message Passing Interface (MPI) Library. 
        However, in doing so, you must properly fence the MPI section with the proper preprocessor macro, as explained in the previous section.  
    +   The minimal and efficient use of preprocessor macros and directives is allowed and encouraged to avoid code duplication and improve productivity.  

+   All Coarray variables or objects (that are communicated between processes) must begin with the `co_` prefix. Example: `co_logFunc`.  
    **Why?** This helps quickly identify objects responsible for inter-process communications in the Coarray implementation of the library.  

+   All **m**odule **v**ariables or module objects must begin with `mv_`, unless there is a reasonable justification for not doing so. Example: `mv_state`, `comv_logFuncState`.  
    Why? Module variables are global, at least within the module's scope and, at most, everywhere in the program. 
    Global variables are susceptible to creating difficult, hard-to-identify bugs in the code and are hard to follow. 
    Adding the prefix `mv_` to all global variables resolves at least part of the problem by explicitly declaring the variable as global to the developer via its name.  

+   All **m**odule **c**onstant variables that should be defined only once throughout the entire runtime must be prefixed with `mc_` instead of `mv_`. Example: `mc_imageCount = 3`.  
    **Why?** Adding the prefix `mc_` helps other developers realize that the object should not change value at runtime beyond the first initialization. 
    Note that this variable type differs from a compile-time constant since its value must be defined at runtime, which could differ from one run to another. 
    However, once set, the value will remain constant throughout the simulation.  

+   All constants (parameters) are upper-case, separated by underscore. Example: `FILE_EXT = ".txt"`.  

+   All module names must begin with the prefix `pm_`, followed by the module name, preferably in `camelCase` style. Example: `pm_kind`.  
    **Why?** The suffix `pm_` stands for the `ParaMonte Module`.  
    Adding this suffix helps create multiple unique entity names from a single base name.  

+   Each module must be stored in a separate source file, and the source file's name must be the same as the module name unless there is an excellent justification to do otherwise. 
    This is a strict rule to maintain the integrity of the library source file and the build scripts.  

+   Submodule must be preferably reserved only for procedure implementation details and named `routines`. Example: `submodule (pm_test) routines`.  

+   If a module has submodules, the name of the submodule file should follow the module name after `@`. For example the following submodule,  
    ```fortran  
    submodule (pm_io) routines
    ...
    end submodule routines
    ```  
    will have to be stored in a source file with the name `pm_io@routines.F90`.  

+   All type and class names should preferably begin with a lowercase letter and, more importantly, must end with the prefix `_type`. Example: `paramcmc_type`.  
    **Why?** This helps create multiple different entity names from a single base name. For example,  
    ```fortran  
    use pm_container, only: css_type
    type(css_type)  :: css
    ```  

#### Case-Sensitivity [⛓⛓⛓⛓](#case-sensitivity-)

+   Although the Fortran programming language is case-INsensitive, please preserve and ensure the case-sensitivity 
    compliance of all object names including module, type, procedure, and variable names in all codes and files.<br>
    This is crucial for the proper hyper-linking of different HTML pages of the ParaMonte library documentation.<br>
    For example, if the module name is `pm_kind`, do not write it in any other way like `PM_kind` or `PM_KIND`.<br>

+   If a variable or procedure argument in a generic interface can be both scalar and array, it should be
    preferably capitalized since a scalar can be envisioned as particular case of an array with one element.

#### Data Hiding and Privacy [⛓⛓⛓⛓](#data-hiding-and-privacy-)

+   Module entities should be all preferably `public` unless an entity is guaranteed to never be needed outside the module,
    in which case, it should be given the `private` attribute. An example of such guaranteed privacy is a type-bound procedure.
    The reason for keeping all module entities `public` is simple. Generic interfaces cannot be passed as dummy arguments.
    Therefore, all module procedures, even though accessible via generic interfaces, should be kept public.

+   Avoid declaring and using global module variables, unless the variable is guaranteed to remain unchanged after program initialization.<br>
    If a module variable is only initialized once throughout the entire life of the program, then it should be preferably prefixed with `mc_`
    to indicate the variable is a <b>m</b>odule runtime <b>c</b>onstant (not a `parameter`, since parameters are compile-time constants).<br>
    in which case, it should be given the `private` attribute. An example of such guaranteed privacy is a type-bound procedure.<br>
    The reason for keeping all module entities `public` is simple. Generic interfaces cannot be passed as dummy arguments.<br>
    Therefore, all module procedures, even though accessible via generic interfaces, should be kept public.<br>

### Controlling the Runtime Checks [⛓⛓⛓](#controlling-the-runtime-checks-)

+   The ParaMonte library has two different checking modes for any build type.<br>
    The build mode is determined (by the choice of user) via the preprocessor macro `CHECK_ENABLED`.<br>
    If this macro is defined, then certain checks on procedure arguments and validity tests in various
    modules and procedures of the library are automatically performed at runtime to ensure the validity
    and accuracy of the calculations. This will however, lower the runtime performance of the library.<br>
    As such, the ParaMonte build scripts define the preprocessor macro `CHECK_ENABLED` by default for
    the `debug` and `testing` build types, and undefine the macro for the `release` build.<br>
    For any build type, this behavior can be overridden by the user/developer when building the library.<br>
    Refer to the help file of the build scripts for more details and instructions.<br>
    >   **WARNING**  
    >   Do **NOT** lower the case of any uppercase `PURE` procedure attribute.<br>
    >   The uppercase `PURE` is a macro that expands to `pure` when runtime checks are disabled,
    >   that is, when `CHECK_ENABLED` is undefined or set to `0`, otherwise, it is empty.<br>
    >   The worst consequence of failing to pay attention to this matter is that the library will not 
    >   compile and the developer will have to fix any mistakes they have made with `PURE` preprocessor macros.

+   When the library is built with tests enabled, a number of fatal `error stop` statements are disabled in 
    some procedures (like the ParaMonte samplers), in particular, when the library is built in parallel mode.<br>
    Such temporary disabling of global program halts is needed to ensure the error signals are propagated correctly 
    with the library and across multiple processors. However, such extra testing communications (in particular, in parallel mode) 
    can be very costly and therefore, **the library should never be built for production release with any tests of the library enabled**.

## Final Steps [⛓⛓](#final-steps-)

Once you have implemented your contributions,  

+   Do not forget to test your contributions by adding new tests to the library's unit-testing framework.  
+   Also, generate a code coverage report to ensure your contributions do not lower the overall code coverage of the library routines.  
+   Follow the [generic contribution guidelines](../../CONTRIBUTING.md/#all-contributors) to submit and merge your contributions with the library's main branch on GitHub. 

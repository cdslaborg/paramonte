/// \mainpage
///
/// \tableofcontents
///
/// This is the <b>ParaMonte C++</b> documentation website for the C++ users and developers.
///
/// [⛓](#ParaMonte)
/// <!--==============================-->
/// \section ParaMonte What is ParaMonte?
/// <!--==============================-->
///
/// ParaMonte is a multi-language library of serial and parallel Monte Carlo and Machine Learning routines
/// scientific inference, e.g., for sampling mathematical density functions of arbitrary-dimensions,
/// with the design goal of unifying
///
/// +   the **automation** of simulations and inference,
/// +   the **user-friendliness** of the library and routines,
/// +   the **accessibility** from multiple programming environments,
/// +   the **high-performance** at runtime, and,
/// +   the **scalability** across many parallel processors.
///
/// [⛓](#ParaMonteRepository)
/// <!--==================================================-->
/// \section ParaMonteRepository ParaMonte Project Repository
/// <!--==================================================-->
///
/// The ParaMonte library is open-source and is permanently located and maintained on **GitHub** at:
///
/// &emsp;<a href="https://github.com/cdslaborg/paramonte" target="_blank"><b>https://github.com/cdslaborg/paramonte</b></a>
///
/// [⛓](#ParaMonteReleases)
/// <!--===============================================-->
/// \section ParaMonteReleases ParaMonte Prebuilt Releases
/// <!--===============================================-->
///
/// The pre-built releases of the ParaMonte library for select configurations and compilers are available on **GitHub Release** page at:
///
/// &emsp;<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><b>https://github.com/cdslaborg/paramonte/releases</b></a>
///
/// For instructions to build the ParaMonte library from source files, visit the ParaMonte library main documentation website [linked below](@ref ParaMonteDocumentation).
///
/// [⛓](#ParaMonteDocumentation)
/// <!--========================================================-->
/// \section ParaMonteDocumentation ParaMonte Documentation Website
/// <!--========================================================-->
///
/// For information about the ParaMonte library in general and in other supported programming languages, visit:
///
/// &emsp;<a href="https://www.cdslab.org/paramonte" target="_blank"><b>https://www.cdslab.org/paramonte</b></a>
///
/// [⛓](#ParaMonteLangDocumentation)
/// \section ParaMonteLangDocumentation ParaMonte C++ Documentation Website
///
/// The documentation for the latest version of the ParaMonte C++ library is always <a href="../../cpp/latest/html/index.html" target="_blank"><b>available on this page</b></a>.<br>
///
/// [⛓](#ParaMonteLangModules)
/// \section ParaMonteLangModules ParaMonte C++ Library Modules
///
/// The ParaMonte C++ library currently contains,
///
/// <ol>
///     <li>    A set of routines that can be used for **optimization and sampling of mathematical density functions with arbitrary precision**.
/// </ol>
///
/// The ParaMonte C++ library is currently under active development to extend the library functionalities to other tasks supported by the [ParaMonte Fortran library](../../../fortran/latest/html).<br>
/// For a full list of all available functionalities and modules, see the [modules listing](./namespaces.html) of this documentation website.<br>
///
/// Module                  | Functionality
/// ------------------------|--------------
/// [pm_sampling](@ref pm_sampling) | This module contains procedures and generic interfaces for the ParaMonte C++ library sampler routines.
///
/// The following modules are currently under development for future addition to the ParaMonte C++ library.
///
/// Module                  | Functionality
/// ------------------------|--------------
/// pm_optimization         | This module contains procedures, generic interfaces, and types for numerical optimizations of mathematical functions.
/// pm_quadPack             | This module contains procedures for non-adaptive and adaptive global numerical quadrature.
/// pm_quadRomb             | This module contains procedures for numerical integrations.
///
/// [⛓](#ParaMonteLangNamingConventions)
/// \section ParaMonteLangNamingConventions ParaMonte C++ Naming Conventions
///
/// +   The **CamelCase** naming style is enforced throughout the ParaMonte C++ library.
///
/// [⛓⛓](#ParaMonteLangNamingConventionsVariables)
/// \subsection ParaMonteLangNamingConventionsVariables ParaMonte C++ Naming Conventions: Variables
///
/// By convention in this library,
///
/// +   All variables and procedure names except compile-time constants begin with an lowercase letter.
///
/// +   All constants and parameters should generally be typed in uppercase.
///
/// +   The names of variables that **always** represent vectors of values **may** be suffixed with `Vec` or `Vector`.
///
/// +   The names of variables that **always** represent matrices of values **can** be suffixed with `Mat` or `Matrix`.
///
/// +   A significant attempt has been made to end all Boolean variables with a passive verb.<br>
///     This is to ensure that the full variable name virtually forms a proposition.<br>
///     In other words, a Boolean variable name should be an English-language statement that evaluates to either `true` or `false`.<br>
///     For example, `parallelismMpiFinalizeEnabled` is one such proposition.<br>
///     +   Occasionally, names that begin with the verb `is` can also be used to label `_Bool` objects.<br>
///     +   But as a general rule, names that begin with a verb should be reserved for procedures.<br>
///
/// [⛓⛓](#ParaMonteLangNamingConventionsProcedures)
/// \subsection ParaMonteLangNamingConventionsProcedures ParaMonte C++ Naming Conventions: Procedures
///
/// +   Procedure names should be descriptive of the action performed by the procedure. For example,
///
///     +   `getCov` means *generate (a) covariance matrix*.
///     +   `setMatCopy` means *copy the matrix into the specified buffer*.
///
/// +   Procedure names should virtually always begin with a lowercase verb.
///
/// +   **Exceptions** to the above and below rules are allowed when,
///     +   the procedure name is exceptionally famous, or
///     +   it is very inconvenient to prefix the procedure name with a verb.
///
/// +   The names of functions with return values should preferably begin with <b>`get`</b>.<br>
///     The reasoning is simple: Functions in mathematics *generate and obtain* a **new** object instead of changing (resetting) the state of an existing object.<br>
///     For example, the function `getMatSym(mat)` generates a symmetric version of the input matrix and returns it as the function result.<br>
///     **Exceptions** to this naming convention are allowed, for example, when the `get` prefix is inconvenient or when the function returns a `_Bool` result.<br>
///
/// +   Functions that return objects of type `_Bool` should be preferably prefixed with `is` or
///     be named such that the name begins with a verb and reads as a proposition, evaluating to
///     either true or false.
///
/// +   The keyword <b>`get`</b> should be avoided as a prefix for `void` function names since `void` functions do not generate
///     and <i>get</i> a new object as their results, rather they <b>(re)set</b> the state of existing objects passed to them.<br>
/// +   A more appropriate prefix for `void` functions is `set`, i.e., `setReplaced()`.
///
/// [⛓](#ParaMonteLangAbbreviationGuidlines)
/// \section ParaMonteLangAbbreviationGuidlines ParaMonte C++ Abbreviation Guidelines
///
/// The following list of abbreviations is in alphabetical order to enable faster search:
///
/// +   The abbreviation `avg`      stands for **average** (rarely used).
/// +   The abbreviation `cdf`      stands for **Cumulative Distribution Function** in the context of statistics. Example: `getNormCDF()`.
/// +   The abbreviation `cho`      stands for **Cholesky factorization**. Example: `setChoLow()`.
/// +   The abbreviation `chol`     stands for **Cholesky factorization**. Example: `setMatChol()`.
/// +   The abbreviation `cor`      stands for **correlation**. Example: `getCor()`.
/// +   The abbreviation `cov`      stands for **covariance**. Example: `getCov()`.
/// +   The abbreviation `cum`      stands for **cumulative**. Example: `getCumSum()`.
/// +   The abbreviation `coef`     stands for **coefficient**. Example: `corcoef_type()`.
/// +   The abbreviation `def`      stands for **default** in variable names (mostly as a prefix `def_` or suffix `_def`).
/// +   The abbreviation `def`      stands for **definite** (mostly in procedure names dealing with positive-definite matrices)
/// +   The abbreviation `den`      stands for **density**, mostly in the context of statistical procedures and objects. Example: `getLogProbDen()`.
/// +   The abbreviation `det`      stands for **determinant**, mostly in the context of Matrix and linear algebra. Example: `getMatDet()`.
/// +   The abbreviation `dia`      stands for **diagonal**, mostly in the context of matrix algebra, matrix packing, or Cholesky factorization. Example: `dia_type()`.
/// +   The abbreviation `diag`     stands for **diagonal**, mostly as dummy argument in matrix algebra procedures.
/// +   The abbreviation `desc`     stands for **description**, mostly as a dummy argument in tests.
/// +   The abbreviation `diff`     stands for **difference**. Example: `setDisSortedExpDiff()`.
/// +   The abbreviation `dist`     stands for **distance** or **distribution** depending on the context. Example: `DistMulti_type`.
/// +   The abbreviation `eff`      stands for **effective**. Example: `effSamSize`.
/// +   The abbreviation `exp`      stands for **exponential** or **exponentiated**. Example: `setDisSortedExpDiff()`.
/// +   The abbreviation `hell`     stands for **Hellinger** in statistical distance computations. Example: `getDisHellSq()`.
/// +   The abbreviation `herm`     stands for **hermitian** in matrix algebra.
/// +   The abbreviation `ICE`      stands for **Internal Compiler Error**. It typically appears in the bug descriptions tagged via Doxygen command <tt>\\bug</tt>.
/// +   The abbreviation `inv`      stands for **inverse**. Example: `getMatInv()`.
/// +   The abbreviation `ks`       stands for **Kolmogorov-Smirnov** test. Example: `getProbKS()`.
/// +   The abbreviation `lin`      stands for **linear**. Example: `getLinSpace()`.
/// +   The abbreviation `low`      stands for **lower triangle of a matrix** or **lower limits**. Example: `setChoLow()`.
/// +   The abbreviation `mahal`    stands for **Mahalanobis** distance. Example: `getDisMahalSq()`.
/// +   The abbreviation `mat`      stands for **matrix**. Example: `getMatInv()`.
/// +   The abbreviation `multi`    stands for **multivariate** mostly used in the context of statistical distributions. Example: `getMultiNormRand()`.
/// +   The abbreviation `msn`      stands for **Multivariate Skew-Normal** mostly used in the context of the statistical MultiVariate Skew-Normal distribution.
/// +   The abbreviation `mvn`      stands for **MultiVariate Normal** mostly used in the context of the statistical MultiVariate Normal distribution.
/// +   The abbreviation `mvu`      stands for **MultiVariate Uniform** mostly used in the context of the statistical MultiVariate (ellipsoidal) Uniform distribution.
/// +   The abbreviation `norm`     stands for **normal** in the context of statistical distributions or **normalization** factor. Example: `DistMultiNorm_type`.
/// +   The abbreviation `normed`   stands for **normalized** mostly in the context of statistical samples. Example: `NormedSample`.
/// +   The abbreviation `pdf`      stands for **Probability Density Function** in the context of statistics. Example: `getNormLogPDF()`.
/// +   The abbreviation `pos`      stands for **positive**. Example: `getInvPosDefMat()`.
/// +   The abbreviation `prob`     stands for `probability`, mostly in the context of statistical applications. Example: `getLogProb()`.
/// +   The abbreviation `proc`     stands for `procedure`, particularly, when it appears as the suffix `_proc` in `abstract interface` definitions.
/// +   The abbreviation `quan`     stands for **quantile**, mostly in the context of statistics. Example: `getParetoLogQuan()`.
/// +   The abbreviation `rand`     stands for **random**, mostly in the context of statistics. Example: `getUnifRand()`.
/// +   The abbreviation `ref`      stands for **reference**, mostly in the context of testings to represent the reference values for comparison. Example: `mean_ref`.
/// +   The abbreviation `sam`      stands for **sample**, mostly in the context of statistics. Example: `effSamSize`.
/// +   The abbreviation `sq`       stands for **squared**. Example: `getDisMahalSq()`.
/// +   The abbreviation `stat`     stands for **statistics**. Example: `StatDRAM_type`.
/// +   The abbreviation `std`      stands for **standard deviation**. Example: `StdVec`.
/// +   The abbreviation `sym`      stands for **symmetric**.
/// +   The abbreviation `symm`     stands for **symmetric**.
/// +   The abbreviation `udf`      stands for **Unnormalized Density Function** in the context of statistics. Example: `getEggBoxLogUDF()`.
/// +   The abbreviation `uni`      stands for **univariate**, mostly used in the context of statistical distributions. Example: `DistUni_type`.
/// +   The abbreviation `unif`     stands for **uniform**, mostly in the context of the uniform statistical distribution. Example: `getUnifRand()`.
/// +   The abbreviation `upp`      stands for **upper triangle of a matrix** or **upper limits**. Example: `setChoUpp()`.
/// +   The abbreviation `vec`      stands for **vector**. Example: `stdVec`.
///
/// [⛓](#ParaMonteLangDeveloperWarnings)
/// \section ParaMonteLangDeveloperWarnings ParaMonte C++ Developer Guidelines and Warnings
///
/// The ParaMonte C++ library development and guidelines are summarized in
/// <a href="https://github.com/cdslaborg/paramonte/blob/main/src/cpp/CONTRIBUTING.md" target="_blank">CONTRIBUTING.md</a>.
///
/// [⛓](#ParaMonteLangDocumentationGuidelines)
/// \section ParaMonteLangDocumentationGuidelines ParaMonte C++ Documentation Guidelines
///
/// +   **Doxygen custom command orderings**.
///
///     +   The Doxygen tag `\brief` must always be the first line of the documentation of modules, types, and procedures.<br>
///         Example: [pm_sampling](@ref pm_sampling).<br>
///     +   The Doxygen tag `\details`, if it exists, must always immediately follow the Doxygen tag `\brief`.<br>
///         Example: [pm_sampling](@ref pm_sampling).<br>
///     +   The Doxygen tag `\param`, if any number of it exists, must always immediately follow the Doxygen tag `\brief` (or `\details` if it exists).<br>
///         Example: [runParaDRAMD()](@ref pm_sampling::runParaDRAMD).<br>
///     +   The Doxygen tag `\return`, must be exclusively used to indicate the return value of functions.<br>
///         If it exists, it must appear immediately after the set of `\param` tags. Example: [runParaDRAMD()](@ref pm_sampling::runParaDRAMD).<br>
///     +   If a generic interface is being documented, the ParaMonte custom command <tt>\\interface</tt> must appear immediately
///         after the Doxygen `\return`, `\param`, `\details`, or `\brief` tags in the specified order, if any exists.<br>
///     +   The Doxygen tag `\warning`, if any number of it exists, must immediately follow the Doxygen tag `\return` if it exists,
///         otherwise `\param` if it exists, otherwise `\details` if it exists, otherwise `\brief`.<br>
///         The `\warning` tag must be used to highlight situations that require special attention of the user,
///         otherwise, there is a danger for the code section being documented to not behave normally as one may expect.<br>
///     +   The Doxygen tag `\attention` has the same functionality and usage as `\warning`.<br>
///         Therefore, `\warning` should be preferred wherever `\attention` is needed.<br>
///         Exceptions are allowed and if they occur, the same documentation conventions as those of `\warning` also apply to the tag `\attention`.<br>
///     +   The Doxygen tag `\remark`, if any number of it exists, must immediately follow the Doxygen tag `\warning` if it exists,
///         otherwise the Doxygen tag `\return` if it exists, otherwise `\param` if it exists, otherwise `\details` if it exists, otherwise `\brief`.<br>
///         The tag `\remark` should be reserved for explaining behavior that is **directly** related to the code segment being documented,
///         but its knowledge is not so critical as warrant the use of a `\warning` tag.<br>
///     +   The Doxygen tag `\note`, if it exists, must appear after all `\warning` and `\attention` and `\remark` tags and
///         immediately before the ParaMonte custom command tag `\see` if it exists, otherwise immediately before <tt>\\example</tt> for examples (if it exists).<br>
///     +   The Doxygen tag `\see`, if it exists, must appear after all `\warning` and `\remark` and `\note` tags.<br>
///         If more than one item for the `\see` command exists, each must be written on a separate line and each line must end with the HTML line-break tag `<br>`.
///         Example: See [below](#example-ParaMonteLangDocumentationGuidelines).<br>
///     +   If any example exists, it must appear immediately after the `\see` tag, otherwise after `\note`, `\remark`, `\warning`, `\param`, `\details`, or `\brief` if any exists.<br>
///         ParaMonte examples are initiated by the custom command <tt>\\example</tt> devised in the `config.txt` file of ParaMonte Doxygen documentation.<br>
///         If the example exists in an external file, then it must be included via the Doxygen `\include` command, followed immediately by
///         the ParaMonte custom Doxygen command <tt>\\compile</tt> which inserts the generic example compile commands for the example, followed optionally
///         but immediately by the output file of the example inserted in the documentation via the `\include` command, followed immediately by the inclusion
///         of any other visualization or postprocessing scripts and output.<br>
///         <b>In all steps, it is imperative to not leave any empty lines between the successive commands of the example
///         section</b>, designated by the <tt>\\example</tt>, otherwise, each empty line will start a new paragraph in the documentation.<br>
///         Example: See [below](#example-ParaMonteLangDocumentationGuidelines).<br>
///     +   The Doxygen `\test` tag, if any exists, must appear immediately after the example section designated by the <tt>\\example</tt> tag.<br>
///     +   The Doxygen `\bug` tag, if any exists, must appear immediately after the `\test` tag or any other tag immediately preceding it.<br>
///     +   The Doxygen `\todo` tag, if any exists, must appear immediately after the `\todo` tag or any other tag immediately preceding it.<br>
///     +   The closing command of each documentation section must be the ParaMonte custom command <tt>\\final</tt> separated from the tags before and after by an empty line.<br>
///     +   The Doxygen `\author` tag is the last command to appear in any documentation section, and it must preferably have the format exemplified in the example below.<br>
///     <br>
///     <br>
///
/// +   **ParaMonte Doxygen custom commands**.<br>
///     To simplify documentation and avoid retyping certain frequently used keywords and sentences, a number of Doxygen aliases are predfined
///     in the ParaMonte Doxygen `config.txt` file. These include (but are not limited to):
///         +   <tt>\\warnpure</tt> Inserts a `\warning` about procedures that are `impure` when the library is built the preprocessor macro `CHECK_ENABLED=1`.
///         +   <tt>\\elemental</tt> Inserts a `\remark` tag indicating that the procedure of interest is `elemental`.
///         +   <tt>\\pure</tt> Inserts a `\remark` tag indicating that the procedure of interest is `pure`.
///         +   <tt>\\interface</tt> Starts a **Possible calling interfaces** paragraph where different calling interfaces of a procedure can be listed.
///         +   <tt>\\benchmark</tt> Starts a new **Benchmark** paragraph which is hyper-linked to the generic anchor `#benchmark` at the same location on the same page.
///         +   <tt>\\benchmark{xxx}</tt> Starts a new **Benchmark** paragraph which is hyper-linked to the specific anchor `#benchmark-xxx` at the same location on the same page.
///         +   <tt>\\benchmark{xxx, This is the benchmark title}</tt> Starts a new **Benchmark** paragraph which is hyper-linked to
///             the specific anchor `#benchmark-xxx` at the same location on the same page with the title `This is the benchmark title`.
///         +   <tt>\\example</tt> Starts a new **Example usage** paragraph which is hyper-linked to the generic anchor `#example` at the same location on the same page.
///         +   <tt>\\example{xxx}</tt> Starts a new **Example usage** paragraph which is hyper-linked to the specific anchor `#example-xxx` at the same location on the same page.
///         +   <tt>\\compile</tt> Inserts the set of example compile commands.
///         +   <tt>\\output</tt> Inserts a title line for the output section of an example paragraph.
///         +   <tt>\\postproc</tt> Inserts a title line for the postprocessing section of an example paragraph.
///         +   <tt>\\abbr</tt> Inserts a `\remark` tag about the naming abbreviations used in the library.
///         +   <tt>\\naming</tt> Inserts a `\remark` tag about the naming conventions used in the library.
///         +   <tt>\\license</tt>  Inserts a `\remark` tag about the generic licensing of the library.
///         +   <tt>\\final</tt>    Inserts the set of final generic remarks that should appear at the end of each documentation section.
///     <br>
///     <br>
///
///     For an up-to-date list of all available aliases, check the value of the Doxygen `ALIASES` option in `config.txt` in the ParaMonte C++ documentation repository.
///
/// +   **Escaping the Doxygen reserved characters**.<br>
///     Doxygen has a set of reserved characters whose usage in the documentation must be handled properly.<br>
///         +   Most importantly, the backslash character `\` begins a Doxygen command.<br>
///             To print a backslash character to the output one should escape it via `\\`.<br>
///         +   Also, the use of the percentage symbol `%` requires special care in some instances.<br>
///             This is particularly important when defining Windows environment variables that should typically be enclosed with percentage character.<br>
///     <br>
///
///     For more information, see [the relevant page on Doxygen documentation website](https://www.doxygen.nl/manual/commands.html#cmdfdollar).
///     <br>
///
/// +   Avoid the insertion of an empty documentation line between any two lines of a single Doxygen paragraph.<br>
///     This is crucial when the whole paragraph is indented by a vertical line as is done by Doxygen for `\warning`, `\remark`, `\note` and other similar tags.<br>
///     \example{ParaMonteLangDocumentationGuidelines}
///     The following is an example documentation for a procedure:
///     \verbatim
///
///
/// \ingroup pm_sampling
///
/// \defgroup runParaDRAM runParaDRAM
///
/// \brief
/// Generate and return a non-zero value (`1`) if the procedure fails to fully accomplish the task of generating
/// a Monte Carlo sample of the specified input mathematical objective function, otherwise, return `0`.
///
/// \details
/// This interface group is the entry point to all **C-style interfaces** to the ParaDRAM samplers of mathematical density functions.<br>
/// Although the procedures of this generic interface return a single scalar of type `int32_t`, the procedures generate
/// massive amounts of information about each simulation which are stored in appropriate external hard drive files.<br>
///
/// See,
/// <ol>
///     <li>    [this generic documentation page](\pmdoc_usage_sampling/paradram/output/)
///             for more information on the generated output files for samplings performed using the [ParaDRAM](@ref runParaDRAM) sampler.
/// </ol>
///
/// \param[in]  getLogFunc  :   The input user-specified procedure pointer to the natural logarithm of the target density function.<br>
///                             On input, the function must take an input vector `state` of size `ndim` of type floating-point of kind \C_RKALL
///                             representing a state (point) from within the domain of the user-specified target density function whose function value must be returned.<br>
///                             On output, the user-specified procedure `getLogFunc()` must return the function value corresponding to the input `state[ndim]`.<br>
///                             The following illustrate the generic interface of input function pointer `getLogFunc(state, ndim)`,
///                             \code{.F90}
///                                 REAL getLogFunc(REAL[] state, int32_t ndim);
///                             \endcode
///                             where `REAL` can be a floating-point type of kind \C_RKALL for the corresponding varying-precision sampler interfaces:<br>
///                             <ol>
///                                 <li>    [runParaDRAMF](@ref runParaDRAMF),
///                                 <li>    [runParaDRAMD](@ref runParaDRAMD),
///                                 <li>    [runParaDRAML](@ref runParaDRAML).
///                             </ol>
/// \param[in]  ndim        :   The input scalar constant of type `int32_t` representing the number of dimensions of the domain of the objective function.
/// \param[in]  input       :   The input scalar pointer of type `char` representing the null-terminated C string
///                             that contains either,
///                             <ol>
///                                 <li>    the path to an external input file containing the namelist group of ParaDRAM sampler specifications
///                                         as outlined in the corresponding page of [ParaMonte library generic documentation website](\pmdoc_usage_sampling/paradram/specifications/).<br>
///                                 <li>    the namelist group of ParaDRAM sampler specifications as the can appear in an external input specification file.<br>
///                             </ol>
///                             While all input simulation specifications are optional, it is highly recommended to pay
///                             attention to the default settings of the domain boundaries and sampler starting point.<br>
///                             (**optional**. It is considered as missing if set to `NULL`.)
///
/// \return
/// `stat`                  :   The output scalar of type `int32_t` that is `0` **if and only if**
///                             the sampler succeeds in sampling the specified density function.<br>
///
/// \interface{runParaDRAM}
/// \code{.cpp}
///
///     include "pm_sampling.hpp"
///
///     stat = runParaDRAMF(getLogFunc, ndim, input)
///     stat = runParaDRAMD(getLogFunc, ndim, input)
///     stat = runParaDRAML(getLogFunc, ndim, input)
///
/// \endcode
///
/// \warning
/// Beware that the definition of extended precision real type `long double` is compiler and platform dependent
/// making the use of `long double` with precompiled ParaMonte libraries problematic and non-functional.<br>
///
/// \warning
/// The condition `0 < ndim` must hold for the corresponding input arguments.<br>
/// \vericon
///
/// \example{mvn}
/// \include{lineno} example/pm_sampling/mvn/main.cpp
/// \inputfile{mvn, input.nml, example specifications namelist input file}
/// \include{lineno} example/pm_sampling/mvn/input.nml
/// \cmakefile{mvn}
/// \include{lineno} example/pm_sampling/mvn/CMakeLists.txt
/// \cmakerun{mvn}
/// \postproc{mvn}
/// \include{lineno} example/pm_sampling/mvn/main.py
/// \vis{mvn}
/// \image html example/pm_sampling/mvn/mvn.traceplot.png width=700
/// \image html example/pm_sampling/mvn/mvn.scatterplot.png width=700
/// \image html example/pm_sampling/mvn/mvn.proposalAdaptation.png width=700
///
/// \example{himmelblau}
/// \include{lineno} example/pm_sampling/himmelblau/main.cpp
/// \inputfile{himmelblau, input.nml, example specifications namelist input file}
/// \include{lineno} example/pm_sampling/himmelblau/input.nml
/// \cmakefile{himmelblau}
/// \include{lineno} example/pm_sampling/himmelblau/CMakeLists.txt
/// \cmakerun{himmelblau}
/// \postproc{himmelblau}
/// \include{lineno} example/pm_sampling/himmelblau/main.py
/// \vis{himmelblau}
/// \image html example/pm_sampling/himmelblau/himmelblau.traceplot.png width=700
/// \image html example/pm_sampling/himmelblau/himmelblau.scatterplot.png width=700
/// \image html example/pm_sampling/himmelblau/himmelblau.proposalAdaptation.png width=700
///
/// \test
/// [test_pm_sampling](@ref test_pm_sampling)
///
/// \bug
/// \status \unresolved
/// \source \ifx{2024.0.0 20231017}
/// \desc
/// \ifx This interface yields a segmentation fault error for all `real` types supported when linked with the ParaMonte library built with \ifx.
/// \remedy
/// For now, only \ifort will be used.<br>
///
/// \bug
/// \status \unresolved
/// \source \ifort{2021.11.0 20231010}
/// \desc
/// The `runParaDRAML` interface for `long double` yields a segmentation fault error.
/// \remedy
/// None as of today.<br>
///
/// \todo
/// \pvhigh
/// The current tests for `long double` interface fail with \ifx and \icx compilers.<br>
/// Apparently, there is a mismatch between `long double` and `c_long_double` storage.<br>
/// This issue does not exist with GNU compilers, although the GNU definition of `long double`
/// appears to yield incorrect values in some calculations, e.g., in `isFailedGeomCyclicFit()` of the ParaMonte library.<br>
///
/// \todo
/// \pmed
/// This C-style interface should be ideally wrapped in a C++ template.
///
/// \final{runParaDRAM}
///
/// \author
/// \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
///
///     \endverbatim
///     <br>
///     The above example documentation snippet will generate [an HTML similar to this documentation](@ref runParaDRAM).<br>
///     Note the lack of an empty line among the commands that immediately follow <tt>\\example</tt>.<br>
///     This is essential to keep the entire example section in the same paragraph.
///
/// <br>
///
/// [⛓](#ParaMonteLangExamples)
/// \section ParaMonteLangExamples ParaMonte C++ Language Examples
///
/// The ParaMonte C++ library ships with tens of thousands of example usage that are available in the `example/cpp` folder in the [root directory of the project repository](https://github.com/cdslaborg/paramonte).<br>
/// These examples are also available and discussed in the documentations of individual modules and procedures of this this documentation website.<br>
///
/// [⛓](#ParaMonteLangBenchmarks)
/// \section ParaMonteLangBenchmarks ParaMonte C++ Language Benchmarks
///
/// The ParaMonte C++ library currently does not ship with any benchmarks.<br>
/// If you would like to see a relevant benchmark currently not included, [discuss it here](https://github.com/cdslaborg/paramonte/discussions)
/// or [raise an issue here](https://github.com/cdslaborg/paramonte/issues) for consideration or volunteer to implement it!<br>
///
/// [⛓](#ParaMonteLangDocumentationTroubleshooting)
/// \section ParaMonteLangDocumentationTroubleshooting ParaMonte C++ Documentation Troubleshooting
///
/// <ol>
///     <li>    **Side navigation pane disappears in some documentation pages.**<br>
///             This issue most likely originates from the interference of browser addons with the documentation.<br>
///             This issue is mostly observed on Firefox browsers.<br>
///             If it occurs, open the page in *browser private mode* or use other (e.g., chrome-based) browsers.<br>
/// </ol>
///
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

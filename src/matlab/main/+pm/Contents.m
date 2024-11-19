%>  \mainpage
%>
%>  \tableofcontents
%>
%>  This is the <b>ParaMonte MATLAB</b> documentation website for the MATLAB users and developers.
%>
%>  [⛓](#ParaMonte)
%>  <!--==============================-->
%>  \section ParaMonte What is ParaMonte?
%>  <!--==============================-->
%>
%>  ParaMonte is a multi-language library of serial and parallel Monte Carlo and Machine Learning routines
%>  scientific inference, e.g., for sampling mathematical density functions of arbitrary-dimensions,
%>  with the design goal of unifying
%>
%>  +   the **automation** of simulations and inference,
%>  +   the **user-friendliness** of the library and routines,
%>  +   the **accessibility** from multiple programming environments,
%>  +   the **high-performance** at runtime, and,
%>  +   the **scalability** across many parallel processors.
%>
%>  [⛓](#ParaMonteRepository)
%>  <!--==================================================-->
%>  \section ParaMonteRepository ParaMonte Project Repository
%>  <!--==================================================-->
%>
%>  The ParaMonte library is open-source and is permanently located and maintained on **GitHub** at:
%>
%>  &emsp;<a href="https://github.com/cdslaborg/paramonte" target="_blank"><b>https://github.com/cdslaborg/paramonte</b></a>
%>
%>  [⛓](#ParaMonteReleases)
%>  <!--===============================================-->
%>  \section ParaMonteReleases ParaMonte Prebuilt Releases
%>  <!--===============================================-->
%>
%>  The pre-built releases of the ParaMonte library for select configurations and compilers are available on **GitHub Release** page at:
%>
%>  &emsp;<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><b>https://github.com/cdslaborg/paramonte/releases</b></a>
%>
%>  For instructions to build the ParaMonte library from source files, visit the ParaMonte library main documentation website [linked below](@ref ParaMonteDocumentation).
%>
%>  [⛓](#ParaMonteDocumentation)
%>  <!--========================================================-->
%>  \section ParaMonteDocumentation ParaMonte Documentation Website
%>  <!--========================================================-->
%>
%>  For information about the ParaMonte library in general and in other supported programming languages, visit:
%>
%>  &emsp;<a href="https://www.cdslab.org/paramonte" target="_blank"><b>https://www.cdslab.org/paramonte</b></a>
%>
%>  [⛓](#ParaMonteLangDocumentation)
%>  \section ParaMonteLangDocumentation ParaMonte MATLAB Documentation Website
%>
%>  The documentation for the latest version of the ParaMonte MATLAB library is always <a href="../../matlab/latest/index.html" target="_blank"><b>available on this page</b></a>.<br>
%>
%>  [⛓](#ParaMonteUsage)
%>  <!--============================================-->
%>  \section ParaMonteUsage ParaMonte MATLAB QuickStart
%>  <!--============================================-->
%>
%>  For more information on the installation, general usage, and examples, visit:
%>
%>  &emsp;<a href="https://www.cdslab.org/paramonte/generic/latest/installation/matlab/" target="_blank"><b>https://www.cdslab.org/paramonte/generic/latest/installation/matlab/</b></a>
%>
%>  To get started with the library once you download it, simply add the path to the ParaMonte library `+pm`
%>  package to access the routines and functionalities available in the library.<br>
%>
%>  \code{.m}
%>
%>      pmlibRootDir = './'; % Change this to the directory containing the ParaMonte MATLAB package `+pm`.
%>      addpath(pmlibRootDir);
%>      pm.lib.version();
%>      doc pm
%>
%>  \endcode
%>
%>  [⛓](#ParaMonteLangModules)
%>  \section ParaMonteLangModules ParaMonte MATLAB Library Modules
%>
%>  The MATLAB equivalent of modules is called **package**.
%>  Similar to modules in other languages, MATLAB packages can be imported to the MATLAB environment
%>  or simply be used as namespace. For example,
%>
%>  \code{.m}
%>      pm.lib.version()
%>  \endcode
%>
%>  The ParaMonte MATLAB library currently contains a myriad of MATLAB *packages*.<br>
%>  The ParaMonte MATLAB library is currently under active development to extend the library functionalities to other tasks supported by the [ParaMonte Fortran library](../../../fortran/latest/index.html).<br>
%>  For a full list of all available functionalities and modules, see the [files listing](./files.html) and [class listing](./annotated.html) of this documentation website.<br>
%>
%>  [⛓⛓](#ParaMonteLangModulesSamplers)
%>  \section ParaMonteLangModulesSamplers ParaMonte MATLAB Library Samplers
%>
%>  Perhaps the most relevant ParaMonte MATLAB modules to the end users are the ParaMonte samplers in the `pm.sampling` package of the library.<br>
%>  See for example, the documentation of the parallel Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo sampler [pm.sampling.Paradram](@ref Paradram).<br>
%>
%>  [⛓](#ParaMonteLangNamingConventions)
%>  \section ParaMonteLangNamingConventions ParaMonte MATLAB Naming Conventions
%>
%>  +   The **camelCase** naming style is enforced throughout the ParaMonte MATLAB library for all variable and function names.<br>
%>  +   The **PascalCase** naming style is enforced throughout the ParaMonte MATLAB library for all class names.<br>
%>      In other words, **all class names begin with a capital letter**.<br>
%>
%>  [⛓⛓](#ParaMonteLangNamingConventionsVariables)
%>  \subsection ParaMonteLangNamingConventionsVariables ParaMonte MATLAB Naming Conventions: Variables
%>
%>  By convention in this library,
%>
%>  +   All variables and procedure names except compile-time constants begin with an lowercase letter.
%>
%>  +   All constants and parameters should generally be typed in uppercase.
%>
%>  +   The names of variables that **always** represent vectors of values **may** be suffixed with `Vec` or `Vector`.
%>
%>  +   The names of variables that **always** represent matrices of values **may** be suffixed with `Mat` or `Matrix`.
%>
%>  +   A significant attempt has been made to end all `logical` variables with a passive verb.<br>
%>      This is to ensure that the full variable name virtually forms a proposition.<br>
%>      In other words, a `logical` variable name should be an English-language statement that evaluates to either `true` or `false`.<br>
%>      For example, `parallelismMpiFinalizeEnabled` is one such proposition.<br>
%>      +   Occasionally, names that begin with the verb `is` can also be used to label `logical` objects.<br>
%>      +   But as a general rule, names that begin with a verb should be reserved for procedures.<br>
%>
%>  [⛓⛓](#ParaMonteLangNamingConventionsProcedures)
%>  \subsection ParaMonteLangNamingConventionsProcedures ParaMonte MATLAB Naming Conventions: Procedures
%>
%>  +   Procedure names should be descriptive of the action performed by the procedure. For example,
%>
%>      +   `getCov` means *generate (a) covariance matrix*.
%>      +   `setMatCopy` means *copy the matrix into the specified buffer*.
%>
%>  +   Procedure names should preferably begin with a lowercase verb.
%>
%>  +   **Exceptions** to the above and below rules are allowed when,
%>      +   the procedure name is exceptionally famous, or
%>      +   it is very inconvenient to prefix the procedure name with a verb.
%>
%>  +   The names of functions with return values should preferably begin with <b>`get`</b>.<br>
%>      The reasoning is simple: Functions in mathematics *generate and obtain* a **new** object instead of changing (resetting) the state of an existing object.<br>
%>      For example, the function `getMatSym(mat)` generates a symmetric version of the input matrix and returns it as the function result.<br>
%>      **Exceptions** to this naming convention are allowed, for example, when the `get` prefix is inconvenient or when the function returns a `logical` result.<br>
%>
%>  +   Functions that return objects of type `logical` should be preferably prefixed with `is` or
%>      be named such that the name begins with a verb and reads as a proposition, evaluating to
%>      either `.true.` or `.false.`.
%>
%>  [⛓](#ParaMonteLangAbbreviationGuidlines)
%>  \section ParaMonteLangAbbreviationGuidlines ParaMonte MATLAB Abbreviation Guidelines
%>
%>  The following list of abbreviations is in alphabetical order to enable faster search:
%>
%>  +   The abbreviation `avg`      stands for **average** (rarely used).
%>  +   The abbreviation `cdf`      stands for **Cumulative Distribution Function** in the context of statistics. Example: `getNormCDF()`.
%>  +   The abbreviation `cho`      stands for **Cholesky factorization**. Example: `setChoLow()`.
%>  +   The abbreviation `chol`     stands for **Cholesky factorization**. Example: `setMatChol()`.
%>  +   The abbreviation `cor`      stands for **correlation**. Example: `getCor()`.
%>  +   The abbreviation `cov`      stands for **covariance**. Example: `getCov()`.
%>  +   The abbreviation `cum`      stands for **cumulative**. Example: `getCumSum()`.
%>  +   The abbreviation `coef`     stands for **coefficient**. Example: `corcoef_type()`.
%>  +   The abbreviation `def`      stands for **default** in variable names (mostly as a prefix `def_` or suffix `_def`).
%>  +   The abbreviation `def`      stands for **definite** (mostly in procedure names dealing with positive-definite matrices)
%>  +   The abbreviation `den`      stands for **density**, mostly in the context of statistical procedures and objects. Example: `getLogProbDen()`.
%>  +   The abbreviation `det`      stands for **determinant**, mostly in the context of Matrix and linear algebra. Example: `getMatDet()`.
%>  +   The abbreviation `dia`      stands for **diagonal**, mostly in the context of matrix algebra, matrix packing, or Cholesky factorization. Example: `dia_type()`.
%>  +   The abbreviation `diag`     stands for **diagonal**, mostly as dummy argument in matrix algebra procedures.
%>  +   The abbreviation `desc`     stands for **description**, mostly as a dummy argument in tests.
%>  +   The abbreviation `diff`     stands for **difference**. Example: `setDisSortedExpDiff()`.
%>  +   The abbreviation `dist`     stands for **distance** or **distribution** depending on the context. Example: `DistMulti_type`.
%>  +   The abbreviation `eff`      stands for **effective**. Example: `effSamSize`.
%>  +   The abbreviation `exp`      stands for **exponential** or **exponentiated**. Example: `setDisSortedExpDiff()`.
%>  +   The abbreviation `hell`     stands for **Hellinger** in statistical distance computations. Example: `getDisHellSq()`.
%>  +   The abbreviation `herm`     stands for **hermitian** in matrix algebra.
%>  +   The abbreviation `ice`      stands for **Internal Compiler Error**. It typically appears in the bug descriptions tagged via Doxygen command <tt>\\bug</tt>.
%>  +   The abbreviation `inv`      stands for **inverse**. Example: `getMatInv()`.
%>  +   The abbreviation `ks`       stands for **Kolmogorov-Smirnov** test. Example: `getProbKS()`.
%>  +   The abbreviation `lin`      stands for **linear**. Example: `getLinSpace()`.
%>  +   The abbreviation `low`      stands for **lower triangle of a matrix** or **lower limits**. Example: `setChoLow()`.
%>  +   The abbreviation `mahal`    stands for **Mahalanobis** distance. Example: `getDisMahalSq()`.
%>  +   The abbreviation `mat`      stands for **matrix**. Example: `getMatInv()`.
%>  +   The abbreviation `multi`    stands for **multivariate** mostly used in the context of statistical distributions. Example: `getMultiNormRand()`.
%>  +   The abbreviation `msn`      stands for **Multivariate Skew-Normal** mostly used in the context of the statistical MultiVariate Skew-Normal distribution.
%>  +   The abbreviation `mvn`      stands for **MultiVariate Normal** mostly used in the context of the statistical MultiVariate Normal distribution.
%>  +   The abbreviation `mvu`      stands for **MultiVariate Uniform** mostly used in the context of the statistical MultiVariate (ellipsoidal) Uniform distribution.
%>  +   The abbreviation `norm`     stands for **normal** in the context of statistical distributions or **normalization** factor. Example: `DistMultiNorm_type`.
%>  +   The abbreviation `normed`   stands for **normalized** mostly in the context of statistical samples. Example: `NormedSample`.
%>  +   The abbreviation `pdf`      stands for **Probability Density Function** in the context of statistics. Example: `getNormLogPDF()`.
%>  +   The abbreviation `pos`      stands for **positive**. Example: `getInvPosDefMat()`.
%>  +   The abbreviation `prob`     stands for `probability`, mostly in the context of statistical applications. Example: `getLogProb()`.
%>  +   The abbreviation `proc`     stands for `procedure`, particularly, when it appears as the suffix `_proc` in `abstract interface` definitions.
%>  +   The abbreviation `quan`     stands for **quantile**, mostly in the context of statistics. Example: `getParetoLogQuan()`.
%>  +   The abbreviation `rand`     stands for **random**, mostly in the context of statistics. Example: `getUnifRand()`.
%>  +   The abbreviation `ref`      stands for **reference**, mostly in the context of testings to represent the reference values for comparison. Example: `mean_ref`.
%>  +   The abbreviation `sam`      stands for **sample**, mostly in the context of statistics. Example: `effSamSize`.
%>  +   The abbreviation `sq`       stands for **squared**. Example: `getDisMahalSq()`.
%>  +   The abbreviation `stat`     stands for **statistics**. Example: `StatDRAM_type`.
%>  +   The abbreviation `std`      stands for **standard deviation**. Example: `StdVec`.
%>  +   The abbreviation `sym`      stands for **symmetric**.
%>  +   The abbreviation `symm`     stands for **symmetric**.
%>  +   The abbreviation `udf`      stands for **Unnormalized Density Function** in the context of statistics. Example: `getEggBoxLogUDF()`.
%>  +   The abbreviation `uni`      stands for **univariate**, mostly used in the context of statistical distributions. Example: `DistUni_type`.
%>  +   The abbreviation `unif`     stands for **uniform**, mostly in the context of the uniform statistical distribution. Example: `getUnifRand()`.
%>  +   The abbreviation `upp`      stands for **upper triangle of a matrix** or **upper limits**. Example: `setChoUpp()`.
%>  +   The abbreviation `vec`      stands for **vector**. Example: `stdVec`.
%>
%>  [⛓](#ParaMonteLangDeveloperWarnings)
%>  \section ParaMonteLangDeveloperWarnings ParaMonte MATLAB Developer Guidelines and Warnings
%>
%>  The ParaMonte MATLAB library development and guidelines are summarized in
%>  <a href="https://github.com/cdslaborg/paramonte/blob/main/src/matlab/CONTRIBUTING.md" target="_blank">CONTRIBUTING.md</a>.
%>
%>  [⛓](#ParaMonteLangDocumentationGuidelines)
%>  \section ParaMonteLangDocumentationGuidelines ParaMonte MATLAB Documentation Guidelines
%>
%>  +   **Doxygen custom command orderings**.
%>
%>      +   The Doxygen tag `\brief` must always be the first line of the documentation of modules, types, and procedures.<br>
%>          Example: [pm.sampling.Paradram](@ref Paradram).<br>
%>      +   The Doxygen tag `\details`, if it exists, must always immediately follow the Doxygen tag `\brief`.<br>
%>          Example: [pm.sampling.Paradram](@ref Paradram).<br>
%>      +   The Doxygen tag `\param`, if any number of it exists, must always immediately follow the Doxygen tag `\brief` (or `\details` if it exists).<br>
%>          Example: [pm.sampling.Paradram](@ref Paradram).<br>
%>      +   The Doxygen tag `\return`, must be exclusively used to indicate the return value of functions.<br>
%>          If it exists, it must appear immediately after the set of `\param` tags. Example: [pm.sampling.Paradram](@ref Paradram).<br>
%>      +   If a generic interface is being documented, the ParaMonte custom command <tt>\\interface</tt> must appear immediately
%>          after the Doxygen `\return`, `\param`, `\details`, or `\brief` tags in the specified order, if any exists.<br>
%>      +   The Doxygen tag `\warning`, if any number of it exists, must immediately follow the Doxygen tag `\return` if it exists,
%>          otherwise `\param` if it exists, otherwise `\details` if it exists, otherwise `\brief`.<br>
%>          The `\warning` tag must be used to highlight situations that require special attention of the user,
%>          otherwise, there is a danger for the code section being documented to not behave normally as one may expect.<br>
%>      +   The Doxygen tag `\attention` has the same functionality and usage as `\warning`.<br>
%>          Therefore, `\warning` should be preferred wherever `\attention` is needed.<br>
%>          Exceptions are allowed and if they occur, the same documentation conventions as those of `\warning` also apply to the tag `\attention`.<br>
%>      +   The Doxygen tag `\remark`, if any number of it exists, must immediately follow the Doxygen tag `\warning` if it exists,
%>          otherwise the Doxygen tag `\return` if it exists, otherwise `\param` if it exists, otherwise `\details` if it exists, otherwise `\brief`.<br>
%>          The tag `\remark` should be reserved for explaining behavior that is **directly** related to the code segment being documented,
%>          but its knowledge is not so critical as warrant the use of a `\warning` tag.<br>
%>      +   The Doxygen tag `\note`, if it exists, must appear after all `\warning` and `\attention` and `\remark` tags and
%>          immediately before the ParaMonte custom command tag `\see` if it exists, otherwise immediately before <tt>\\example</tt> for examples (if it exists).<br>
%>      +   The Doxygen tag `\see`, if it exists, must appear after all `\warning` and `\remark` and `\note` tags.<br>
%>          If more than one item for the `\see` command exists, each must be written on a separate line and each line must end with the HTML line-break tag `<br>`.
%>          Example: See [below](#example-ParaMonteLangDocumentationGuidelines).<br>
%>      +   If any example exists, it must appear immediately after the `\see` tag, otherwise after `\note`, `\remark`, `\warning`, `\param`, `\details`, or `\brief` if any exists.<br>
%>          ParaMonte examples are initiated by the custom command <tt>\\example</tt> devised in the `config.txt` file of ParaMonte Doxygen documentation.<br>
%>          If the example exists in an external file, then it must be included via the Doxygen `\include` command, followed immediately by
%>          the ParaMonte custom Doxygen command <tt>\\compile</tt> which inserts the generic example compile commands for the example, followed optionally
%>          but immediately by the output file of the example inserted in the documentation via the `\include` command, followed immediately by the inclusion
%>          of any other visualization or postprocessing scripts and output.<br>
%>          <b>In all steps, it is imperative to not leave any empty lines between the successive commands of the example
%>          section</b>, designated by the <tt>\\example</tt>, otherwise, each empty line will start a new paragraph in the documentation.<br>
%>          Example: See [below](#example-ParaMonteLangDocumentationGuidelines).<br>
%>      +   The Doxygen `\test` tag, if any exists, must appear immediately after the example section designated by the <tt>\\example</tt> tag.<br>
%>      +   The Doxygen `\bug` tag, if any exists, must appear immediately after the `\test` tag or any other tag immediately preceding it.<br>
%>      +   The Doxygen `\todo` tag, if any exists, must appear immediately after the `\todo` tag or any other tag immediately preceding it.<br>
%>      +   The closing command of each documentation section must be the ParaMonte custom command <tt>\\final</tt> separated from the tags before and after by an empty line.<br>
%>      +   The Doxygen `\author` tag is the last command to appear in any documentation section, and it must preferably have the format exemplified in the example below.<br>
%>      <br>
%>      <br>
%>
%>  +   **ParaMonte Doxygen custom commands**.<br>
%>      To simplify documentation and avoid retyping certain frequently used keywords and sentences, a number of Doxygen aliases are predfined
%>      in the ParaMonte Doxygen `config.txt` file. These include (but are not limited to):
%>          +   <tt>\\warnpure</tt> Inserts a `\warning` about procedures that are `impure` when the library is built the preprocessor macro `CHECK_ENABLED=1`.
%>          +   <tt>\\elemental</tt> Inserts a `\remark` tag indicating that the procedure of interest is `elemental`.
%>          +   <tt>\\pure</tt> Inserts a `\remark` tag indicating that the procedure of interest is `pure`.
%>          +   <tt>\\interface</tt> Starts a **Possible calling interfaces** paragraph where different calling interfaces of a procedure can be listed.
%>          +   <tt>\\benchmark</tt> Starts a new **Benchmark** paragraph which is hyper-linked to the generic anchor `#benchmark` at the same location on the same page.
%>          +   <tt>\\benchmark{xxx}</tt> Starts a new **Benchmark** paragraph which is hyper-linked to the specific anchor `#benchmark-xxx` at the same location on the same page.
%>          +   <tt>\\benchmark{xxx, This is the benchmark title}</tt> Starts a new **Benchmark** paragraph which is hyper-linked to
%>              the specific anchor `#benchmark-xxx` at the same location on the same page with the title `This is the benchmark title`.
%>          +   <tt>\\example</tt> Starts a new **Example usage** paragraph which is hyper-linked to the generic anchor `#example` at the same location on the same page.
%>          +   <tt>\\example{xxx}</tt> Starts a new **Example usage** paragraph which is hyper-linked to the specific anchor `#example-xxx` at the same location on the same page.
%>          +   <tt>\\compile</tt> Inserts the set of example compile commands.
%>          +   <tt>\\output</tt> Inserts a title line for the output section of an example paragraph.
%>          +   <tt>\\postproc</tt> Inserts a title line for the postprocessing section of an example paragraph.
%>          +   <tt>\\abbr</tt> Inserts a `\remark` tag about the naming abbreviations used in the library.
%>          +   <tt>\\naming</tt> Inserts a `\remark` tag about the naming conventions used in the library.
%>          +   <tt>\\license</tt> Inserts a `\remark` tag about the generic licensing of the library.
%>          +   <tt>\\final</tt> Inserts the set of final generic remarks that should appear at the end of each documentation section.
%>      <br>
%>      <br>
%>
%>      For an up-to-date list of all available aliases, check the value of the Doxygen `ALIASES` option in `config.txt` in the ParaMonte MATLAB documentation repository.
%>
%>  +   **Escaping the Doxygen reserved characters**.<br>
%>      Doxygen has a set of reserved characters whose usage in the documentation must be handled properly.<br>
%>          +   Most importantly, the backslash character `\` begins a Doxygen command.<br>
%>              To print a backslash character to the output one should escape it via `\\`.<br>
%>          +   Also, the use of the percentage symbol `%` requires special care in some instances.<br>
%>              This is particularly important when defining Windows environment variables that should typically be enclosed with percentage character.<br>
%>      <br>
%>
%>      For more information, see [the relevant page on Doxygen documentation website](https://www.doxygen.nl/manual/commands.html#cmdfdollar).
%>      <br>
%>
%>  +   Avoid the insertion of an empty documentation line between any two lines of a single Doxygen paragraph.<br>
%>      This is crucial when the whole paragraph is indented by a vertical line as is done by Doxygen for `\warning`, `\remark`, `\note` and other similar tags.<br>
%>      \example{ParaMonteLangDocumentationGuidelines}
%>      The following is an example documentation for a procedure:
%>      \verbatim
%>
%>  \brief
%>  This is the ParaDRAM class for generating instances of serial and parallel
%>  Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo
%>  sampler of the ParaMonte MATLAB library.<br>
%>
%>  \brief
%>  Once you assign the desired simulation specifications to the corresponding
%>  attributes within the component `spec` of an object of class [pm.sampling.Paradram](@ref Paradram),
%>  call the ParaDRAM sampler via the object method [pm.sampling.Paradram.run()](@ref Paradram::run).<br>
%>
%>  While the constructor of this class does not take any input arguments,
%>  all ParaDRAM simulation specifications can be set after creating the object.<br>
%>
%>  \return
%>  ``sampler`` :   The output scalar object of class [pm.sampling.Paradram](@ref Paradram).<br>
%>
%>  \interface{Paradram}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>
%>  \endcode
%>
%>  \warning
%>  When using the ParaMonte MATLAB library functionalities, particularly ParaMonte samplers in parallel,
%>  it would be best to close any such aggressive software/applications as Dropbox, ZoneAlarm, ...
%>  that interfere with your ParaMonte MATLAB library output files, potentially causing the tasks
%>  to fail and crash before successful completion.<br>
%>  These situations happen only scarcely.<br>
%>
%>  \note
%>  On Windows systems, when restarting an old interrupted ParaDRAM simulation,
%>  ensure your MATLAB session is also restarted before the simulation restart.<br>
%>  This may be needed as Windows sometimes locks access to some or all of the simulation output files.<br>
%>
%>  \note
%>  To unset an already-set input simulation specification, simply set the
%>  simulation attribute to empty double `[]` or re-instantiate the object.<br>
%>
%>  \see
%>  [ParaDRAM simulation specifications listing](\pmdoc_usage_sampling/paradram/specifications/)<br>
%>  [ParaDRAM simulation restart functionality](\pmdoc_usage_sampling/paradram/restart/)<br>
%>  [ParaDRAM simulation output files](\pmdoc_usage_sampling/paradram/output/)<br>
%>
%>  Example Usage: Serial
%>  ---------------------
%>
%>  First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.
%>
%>  Here is a MATLAB script ``main.m`` for a serial ParaDRAM simulation.<br>
%>  Copy and paste the following code into your MATLAB session:<br>
%>
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      sampler.run ( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function.
%>                  , 4                 ... ndim:       the number of dimensions of the objective function.
%>                  );
%>      samples = sampler.readSample();
%>      sample = samples{1};
%>      tile = pm.vis.TileLine(sample.df);
%>      tile.make("coly", sample.slfc + 1 : sample.slfc + 4, "colc", "sampleLogFunc");
%>
%>  \endcode
%>
%>  The mathematical objective function in the above example is a
%>  is a multivariate Normal distribution centered at the origin,
%>  whose natural logarithm is returned by the lambda (``Anonymous``)
%>  function defined as a function handle input to the ParaDRAM sampler.<br>
%>
%>  Running this code will generate a set of simulation output files (in the current working directory of MATLAB).<br>
%>  Among these, the file suffixed with "_report.txt" contains the full description of all input specifications
%>  of the ParaDRAM simulation as well as other information about the simulation results.<br>
%>
%>  Example Usage: Thread-Parallel
%>  ------------------------------
%>
%>  First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.<br>
%>
%>  Threading parallelism is possible as of ParaMonte MATLAB version ``3.0.0``.<br>
%>  However, only ``singleChain`` ParaDRAM simulations are supported.<br>
%>
%>  Here is a MATLAB script ``main.m`` for a thread-parallel ParaDRAM simulation.<br>
%>  Copy and paste the following code and paste into your MATLAB session:<br>
%>
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      sampler.spec.parallelismNumThread = 0; % use all available threads.
%>      sampler.run ( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function.
%>                  , 4                 ... ndim:       the number of dimensions of the objective function.
%>                  );
%>      samples = sampler.readSample();
%>      sample = samples{1};
%>      pm.vis.tile(sample.contents)
%>
%>  \endcode
%>
%>  The mathematical objective function in the above example is a
%>  is a multivariate Normal distribution centered at the origin,
%>  whose natural logarithm is returned by the lambda (``Anonymous``)
%>  function defined as a function handle input to the ParaDRAM sampler.<br>
%>
%>  Running this code will generate a set of simulation output files (in the current working directory of MATLAB).<br>
%>  Among these, the file suffixed with ``"_report.txt"`` contains the full description of all input specifications
%>  of the ParaDRAM simulation as well as other information about the simulation results.<br>
%>
%>  Specifying ``0`` as the number of threads will lead to using
%>  all available CPU threads for thread-parallel ParaDRAM simulation.<br>
%>
%>  \note
%>  **Benefits of thread-parallelism**<br>
%>  Thread-parallel simulations offer a much more flexible
%>  and easier approach to benefiting from parallelism without
%>  going through the hassle of MPI-parallel simulations.<br>
%>  But they can still potentially offer much faster speed than serial simulations.<br>
%>  The actual speedup depends on a lot of factors.<br>
%>  Moreover, the number of threads is limited to maximum
%>  number of physical cores available on your system.<br>
%>  As such, thread-parallel simulations are not scalable.<br>
%>  If you need scalability, checkout MPI-parallelism below.<br>
%>
%>  Example Usage: MPI-Parallel
%>  ---------------------------
%>
%>  First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.<br>
%>
%>  MPI-parallel simulations can be slightly more cumbersome than thread-parallel simulations
%>  described above because MPI-parallel simulations cannot be performed from within a MATLAB GUI
%>  session and require launching MATLAB via a compatible ``mpiexec`` launcher.<br>
%>
%>  <ol>
%>      <li>    Ensure you need and will get a speedup by running the an MPI-parallel simulation.<br>
%>              Typically, your simulation may then benefit from parallelism only if a single
%>              evaluation of the objective function takes longer than a few milliseconds.<br>
%>
%>      <li>    Ensure the required MPI libraries are installed on your system
%>              (You can skip this step if you know that you already have
%>              a compatible MPI library installed on your system).<br>
%>              On the MATLAB command line type the following,<br>
%>              \code{.m}
%>                  pm.lib.verify();
%>              \endcode
%>              This will verify the existence of a valid MPI library on your system and,
%>              if missing, will guide you to install the MPI library on your system.<br>
%>
%>      <li>    Once the MPI installation is verified, copy and paste the following
%>              code into your MATLAB session:
%>              \code{.m}
%>
%>                  fid = fopen("main_mpi.m", "w");
%>                  sourceCode = ...
%>                  "sampler = pm.sampling.Paradram();" + newline + ...
%>                  "%sampler.mpiname = ''; % set this to an explicit MPI library name if needed." + newline + ...
%>                  "sampler.run( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function." + newline + ...
%>                  "           , 4                 ... ndim:       the number of dimensions of the objective function" + newline + ...
%>                  "           );";
%>                  fprintf(fid, "%s\n", sourceCode);
%>                  fclose(fid);
%>
%>              \endcode
%>
%>      <li>    This will generate a ``main_mpi.m`` MATLAB script file in the current working directory of your MATLAB session.<br>
%>              Now, you can execute this MATLAB script file (``main_mpi.m``) in parallel.<br>
%>              To do so, you need to call MATLAB on a command-line, **out of MATLAB GUI**.<br>
%>              <ol>
%>                  <li>    **On Windows**:<br>
%>                          From within command prompt that recognizes both MATLAB and ``mpiexec``,
%>                          ideally, the Intel dedicated command-prompt that is shipped with Intel MPI library,
%>                          type the following,
%>                          \code{.m}
%>
%>                              mpiexec -localonly -n 3 matlab -batch "main_mpi"
%>
%>                          \endcode
%>
%>                          \note
%>                          In the above MPI launcher command for Windows OS,
%>                          we assumed that you would be using the Intel MPI library, hence,
%>                          the reason for the extra flag ``-localonly``.<br>
%>                          This flag runs the parallel code only on one node, but in doing so,
%>                          it avoids the use of Hydra service and its registration.<br>
%>                          If you are not on a Windows cluster, (e.g., you are using your personal device),
%>                          then we recommend specifying this flag.<br>
%>
%>                  <li>    **On macOS/Linux**:<br>
%>                          From within a Bash terminal that recognizes both MATLAB and ``mpiexec``,
%>                          type the following,
%>                          \code{.m}
%>
%>                              mpiexec -n 3 matlab -batch "main_mpi"
%>
%>                          \endcode
%>
%>                          \note
%>                          In both cases in the above, the script ``main_mpi.m`` will run on 3 processors.<br>
%>                          Feel free to change the number of processors to any number desired.<br>
%>                          But do not request more than the available number of physical cores on your system.<br>
%>              </ol>
%>  </ol>
%>
%>  \warning
%>  Do not add postprocessing codes (such as reading and plotting the output samples) in your MPI-parallel code.<br>
%>  There is no point in doing so, since MATLAB will run in ``-batch`` mode for parallel simulations, disabling all plotting capabilities.<br>
%>  Moreover, if you read and postprocess the output files in parallel mode, the task will be done
%>  by all parallel processes, potentially overwriting external IO activities of each other.<br>
%>  Only perform the sampling as described above in MPI-parallel mode.<br>
%>
%>  ParaDRAM Simulation Specifications
%>  ----------------------------------
%>
%>  The ParaDRAM simulation specifications have lengthy comprehensive descriptions
%>  that appear in full in the output report files of every ParaDRAM simulation.<br>
%>
%>  The best way to learn about individual ParaDRAM simulation attributes
%>  is to a run a minimal serial simulation as given in the above.<br>
%>  You can also use the ``sampler.spec.doc()`` method:
%>
%>  \code{.m}
%>      sampler = pm.sampling.Paradram();
%>      sampler.spec.doc();
%>  \endcode
%>
%>  \final{Paradram}
%>
%>  \author
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>
%>      \endverbatim
%>      <br>
%>      The above example documentation snippet will generate [an HTML similar to this documentation](@ref Paradram).<br>
%>      Note the lack of an empty line among the commands that immediately follow <tt>\\example</tt>.<br>
%>      This is essential to keep the entire example section in the same paragraph.
%>
%>  <br>
%>
%>  [⛓](#ParaMonteLangExamples)
%>  \section ParaMonteLangExamples ParaMonte MATLAB Language Examples
%>
%>  The ParaMonte MATLAB library ships with hundreds of example usage that are available in the `example/matlab` folder in the [root directory of the project repository](https://github.com/cdslaborg/paramonte).<br>
%>  These examples are also available and discussed in the documentations of individual modules and procedures of this this documentation website.<br>
%>
%>  [⛓](#ParaMonteLangBenchmarks)
%>  \section ParaMonteLangBenchmarks ParaMonte MATLAB Language Benchmarks
%>
%>  The ParaMonte MATLAB library currently does not ship with any benchmarks.<br>
%>  If you would like to see a relevant benchmark currently not included, [discuss it here](https://github.com/cdslaborg/paramonte/discussions)
%>  or [raise an issue here](https://github.com/cdslaborg/paramonte/issues) for consideration or volunteer to implement it!<br>
%>
%>  [⛓](#ParaMonteLangDocumentationTroubleshooting)
%>  \section ParaMonteLangDocumentationTroubleshooting ParaMonte MATLAB Documentation Troubleshooting
%>
%>  <ol>
%>      <li>    **Side navigation pane disappears in some documentation pages.**<br>
%>              This issue most likely originates from the interference of browser addons with the documentation.<br>
%>              This issue is mostly observed on Firefox browsers.<br>
%>              If it occurs, open the page in *browser private mode* or use other (e.g., chrome-based) browsers.<br>
%>  </ol>
%>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
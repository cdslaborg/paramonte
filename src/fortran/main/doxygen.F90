!   \mainpage ParaMonte: Parallel Monte Carlo and Machine Learning Library
!>  \mainpage
!>
!>  \tableofcontents
!>
!>  This is the <b>ParaMonte Fortran</b> documentation website for the Fortran users and developers.
!>
!>  [⛓](#ParaMonte)
!>  <!--==============================-->
!>  \section ParaMonte What is ParaMonte?
!>  <!--==============================-->
!>
!>  ParaMonte is a library of serial and parallel Monte Carlo and Machine Learning routines
!>  scientific inference, e.g., for sampling mathematical density functions of arbitrary-dimensions,
!>  with the design goal of unifying
!>
!>  +   the **automation** of simulations and inference,
!>  +   the **user-friendliness** of the library and routines,
!>  +   the **accessibility** from multiple programming environments,
!>  +   the **high-performance** at runtime, and,
!>  +   the **scalability** across many parallel processors.
!>
!>  [⛓](#ParaMonteRepository)
!>  <!--==================================================-->
!>  \section ParaMonteRepository ParaMonte Project Repository
!>  <!--==================================================-->
!>
!>  The ParaMonte library is open-source and is permanently located and maintained on **GitHub** at:
!>
!>  &emsp;<a href="https://github.com/cdslaborg/paramonte" target="_blank"><b>https://github.com/cdslaborg/paramonte</b></a>
!>
!>  [⛓](#ParaMonteReleases)
!>  <!--===============================================-->
!>  \section ParaMonteReleases ParaMonte Prebuilt Releases
!>  <!--===============================================-->
!>
!>  The pre-built releases of the ParaMonte library for select configurations and compilers are available on **GitHub Release** page at:
!>
!>  &emsp;<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><b>https://github.com/cdslaborg/paramonte/releases</b></a>
!>
!>  For instructions to build the ParaMonte library from source files, visit the ParaMonte library main documentation website [linked below](@ref ParaMonteDocumentation).
!>
!>  [⛓](#ParaMonteDocumentation)
!>  <!--========================================================-->
!>  \section ParaMonteDocumentation ParaMonte Documentation Website
!>  <!--========================================================-->
!>
!>  For information about the ParaMonte library in general and in other supported programming languages, visit:<br>
!>
!>  &emsp;<a href="https://www.cdslab.org/paramonte" target="_blank"><b>https://www.cdslab.org/paramonte</b></a>
!>
!>  [⛓](#ParaMonteLangDocumentation)
!>  \section ParaMonteLangDocumentation ParaMonte Fortran Documentation Website
!>
!>  The documentation for the latest version of the ParaMonte Fortran library is always <a href="../../fortran/latest/html/index.html" target="_blank"><b>available on this page</b></a>.<br>
!>
!>  [⛓](#ParaMonteLangModules)
!>  \section ParaMonteLangModules ParaMonte Fortran Library Modules
!>
!>  The ParaMonte Fortran library contains,
!>
!>  <ol>
!>      <li>    a comprehensive list of routines that can be used for **type coercion**.
!>      <li>    a comprehensive list of routines that can be used for **string manipulation**.
!>      <li>    a comprehensive list of routines that can be used for **polynomial calculations**.
!>      <li>    a comprehensive list of routines that can be used for **common filesystem tasks**.
!>      <li>    a comprehensive list of routines that can be used for **calendrical calculations**.
!>      <li>    a comprehensive list of routines that can be used for **common operating system tasks**.
!>      <li>    a comprehensive list of routines that can be used for **computing properties of statistical samples**.
!>      <li>    a comprehensive list of routines that can be used for **common numerical algebraic and mathematical tasks**.
!>      <li>    a comprehensive list of routines that can be used for **manipulating arrays of arbitrary intrinsic types and kind parameters**.
!>      <li>    a comprehensive list of routines that can be used for **computing properties of famous statistical distributions with arbitrary precision**.
!>      <li>    a comprehensive list of routines that can be used for **generating deterministic pseudo-random numbers from various statistical distributions**.
!>      <li>    a comprehensive list of routines that can be used for **optimization and sampling of mathematical density functions with arbitrary precision**.
!>      <li>    a modernization and significant extension of the venerable **QuadPack** Fortran library to support computations with arbitrary precision.
!>      <li>    a modernization and significant extension of the venerable **FFTPACK** Fortran library to support computations with arbitrary precision.
!>      <li>    many other functionalities that are further discussed below or in this documentation.
!>  </ol>
!>
!>  The following is an incomplete list of the functionalities available in the ParaMonte Fortran library.<br>
!>  For a full list of all available functionalities and modules, see the [modules listing](./namespaces.html) of this documentation website.<br>
!>
!>  Module                  | Functionality
!>  ------------------------|--------------
!>  pm_array                | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for various array operations.
!>  pm_arrayCenter          | This module contains procedures and generic interfaces for resizing an input array and centering the original contents of the array in a new array.
!>  pm_arrayChange          | This module contains procedures and generic interfaces for selecting uniformly-distributed random choices from a given character or integer range.
!>  pm_arrayChoice          | This module contains procedures and generic interfaces for selecting uniformly-distributed or arbitrarily-distributed random choices from a given list of intrinsic type of arbitrary kind.
!>  pm_arrayCompact         | This module contains procedures and generic interfaces for condensing (removing duplicate sequential the elements of) an array of arbitrary intrinsic type.
!>  pm_arrayCompareLex      | This module contains procedures and generic interfaces for performing lexicographic comparisons of two arrays of similar type, kind, and rank.
!>  pm_arrayComplement      | This module contains procedures and generic interfaces for computing the absolute or relative complement of one set in another set.
!>  pm_arrayCopy            | This module contains procedures and generic interfaces for copying strided or indexed elements of one scalar string or vector of arbitrary intrinsic type and kind to strided or indexed elements of another scalar string or vector of the same type and kind.
!>  pm_arrayFill            | This module contains procedures and generic interfaces for convenient allocation and filling of arrays of arbitrary intrinsic types (i.e., character, integer, logical, complex, real), kinds, and non-zero ranks (up to 3).
!>  pm_arrayFind            | This module contains procedures and generic interfaces for finding locations of a pattern in arrays of various types at the specified instances of occurrence of pattern.
!>  pm_arrayInit            | This module contains procedures and generic interfaces for efficient initialization of arbitrary rectangular cores and surrounding halos of arrays of arbitrary size, shape, and rank of arbitrary intrinsic type and kind.
!>  pm_arrayInsert          | This module contains procedures and generic interfaces for inserting an insertion into the specified locations of an input arrays of various types.
!>  pm_arrayMembership      | This module contains procedures and generic interfaces for assessing whether particular value(s) or any values or all values within a collection are members of another collection of values, or within a range of values that specifies a mathematical set.
!>  pm_arrayMerge           | This module contains procedures and generic interfaces for sorting and merging two previously-sorted arrays.
!>  pm_arrayMinMax          | This module contains procedures and generic interfaces for finding the minimum and maximum of two input scalar numbers through lexical comparison.
!>  pm_arrayPad             | This module contains procedures and generic interfaces for resizing an input array and padding them with symbols on the left or right.
!>  pm_arrayRange           | This module contains procedures and generic interfaces for generating ranges of discrete character, integer, or real -valued sequences with minimum-possible or user-specified fixed linear spacings.
!>  pm_arrayRank            | This module contains procedures and generic interfaces for obtaining the Ordinal Ranking of the elements of arrays of various types.
!>  pm_arrayRebill          | This module contains procedures and generic interfaces for resizing allocatable arrays of various types, relocating their contents, and rebinding (re-indexing) their lower and upper bounds, and refilling the newly added elements.
!>  pm_arrayRebind          | This module contains procedures and generic interfaces for resizing allocatable arrays of various types, relocating their contents and rebinding (re-indexing) their lower and upper bounds.
!>  pm_arrayRefill          | This module contains procedures and generic interfaces for resizing allocatable arrays of various types, relocating their contents and filling the newly added elements with specific values.
!>  pm_arrayRefine          | This module contains procedures and generic interfaces for refining (thinning) (weighted) arrays of arbitrary intrinsic types.
!>  pm_arrayRemap           | This module contains procedures and generic interfaces for remapping arrays of various types.
!>  pm_arrayRemove          | This module contains procedures and generic interfaces for removing a pattern from arrays of various types at the specified instances of occurrence of pattern.
!>  pm_arrayReplace         | This module contains procedures and generic interfaces for replacing patterns within arrays of various types.
!>  pm_arrayResize          | This module contains procedures and generic interfaces for resizing allocatable arrays of various types and relocating their contents, without initializing or filling the newly added elements with specific values.
!>  pm_arrayReverse         | This module contains procedures and generic interfaces for reversing the order of elements in arrays of various types.
!>  pm_arraySearch          | This module contains procedures and generic interfaces for finding the specific array index whose element has the largest value smaller than the input value in arrays of various types.
!>  pm_arraySelect          | This module contains procedures and generic interfaces for selecting the kth smallest element in unsorted arrays of various types.
!>  pm_arrayShuffle         | This module contains procedures and generic interfaces for shuffling arrays of various types.
!>  pm_arraySort            | This module contains procedures and generic interfaces for various sorting tasks.
!>  pm_arraySpace           | This module contains procedures and generic interfaces for generating arrays with linear or logarithmic spacing.
!>  pm_arraySplit           | This module contains procedures and generic interfaces for splitting arrays of various types at the specified instances of occurrence of pattern.
!>  pm_arrayStrip           | This module contains procedures and generic interfaces for stripping a given pattern from the left and right ends of an array of arbitrary intrinsic type and kind.
!>  pm_arrayUnique          | This module contains procedures and generic interfaces for finding unique values of an input array of various types.
!>  pm_arrayVerbose         | This module contains procedures and generic interfaces for flattening (duplicating the elements of) an array according to a user-specified weight.
!>  pm_batse                | This module contains procedures and generic interfaces for modeling data and detectors of the BATSE Gamma-Ray satellite onboard the NASA Compton Gamma-Ray Observatory.
!>  pm_bench                | This module contains abstract interfaces and types that facilitate benchmarking of different procedures.
!>  pm_bit                  | This module contains constants and procedures that are relevant to bit manipulation.
!>  pm_blas                 | This module contains a set of generic interfaces to the BLAS routines used within the ParaMonte library.
!>  pm_clustering           | This module contains procedures and routines for the computing the Kmeans clustering of a given set of data.
!>  pm_complexAbs           | This module contains procedures and generic interfaces for performing element-wise comparison of the real and imaginary components of scalars and arrays of arbitrary ranks of various types.
!>  pm_complexCompareAll    | This module contains procedures and generic interfaces for checking if both of the corresponding real and imaginary components of two complex numbers satisfy a relational operator.
!>  pm_complexCompareAny    | This module contains procedures and generic interfaces for checking if either of the corresponding real and imaginary components of two complex numbers satisfy a relational operator.
!>  pm_complexCompareLex    | This module contains procedures and generic interfaces for checking if a complex number is lexicographically comparable to another complex number of the same kind.
!>  pm_complexDiv           | This module contains procedures and generic interfaces for computing the complex division robustly without potential overflow of computations.
!>  pm_complexMinMax        | This module contains procedures and generic interfaces for computing element-wise minimum/maximum value/location of the real and imaginary components of scalars and arrays of arbitrary ranks of type complex of arbitrary kinds.
!>  pm_container            | This module contains the derived types for generating allocatable containers of scalar, vector, matrix, or cube of integer, real, complex, logical, and string values of arbitrary kinds.
!>  pm_control              | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.
!>  pm_cosmicRate           | This module contains procedures and generic interfaces for computing the cosmic rates of celestial phenomena.
!>  pm_cosmology            | This module contains procedures and generic interfaces and constants for cosmological calculations.
!>  pm_dateTime             | This module contains classes and procedures for computing, manipulating, and styling dates and times.
!>  pm_distanceBhat         | This module contains classes and procedures for computing the Bhattacharyya statistical distance between two probability distributions.
!>  pm_distanceEuclid       | This module contains procedures and generic interfaces for computing the Euclidean norm of a single point (with respect to origin or a given reference) or the pairwise Euclidean distances (squared) of a collection of points with respect to another set of reference points, optionally without undue overflow or underflow.
!>  pm_distanceHellinger    | This module contains classes and procedures for computing the Hellinger statistical distance between two probability distributions.
!>  pm_distanceKolm         | This module contains classes and procedures for computing the Kolmogorov statistical distance.
!>  pm_distanceMahal        | This module contains classes and procedures for computing the Mahalanobis statistical distance.
!>  pm_distBand             | This module contains procedures and generic interfaces for computing the Band photon distribution widely used in modeling the spectra of a class of celestial objects knowns Gamma-Ray Bursts.
!>  pm_distBern             | This module contains classes and procedures for generating Bernoulli-distributed random numbers.
!>  pm_distBeta             | This module contains classes and procedures for computing various statistical quantities related to the Beta distribution.
!>  pm_distCosRaised        | This module contains classes and procedures for computing various statistical quantities related to the Raised Cosine distribution.
!>  pm_distCov              | This module contains classes and procedures for generating random matrices distributed on the space of positive definite matrices, such that their determinants is uniformly or power-law distributed.
!>  pm_distEggBox           | This module contains classes and procedures for computing various statistical quantities related to the mathematical EggBox density function.
!>  pm_distExp              | This module contains classes and procedures for computing various statistical quantities related to the Exponential distribution.
!>  pm_distExpGamma         | This module contains classes and procedures for computing various statistical quantities related to the ExpGamma distribution.
!>  pm_distGamma            | This module contains classes and procedures for computing various statistical quantities related to the Gamma distribution.
!>  pm_distGenExpGamma      | This module contains classes and procedures for computing various statistical quantities related to the GenExpGamma distribution.
!>  pm_distGenGamma         | This module contains classes and procedures for computing various statistical quantities related to the GenGamma distribution.
!>  pm_distGeom             | This module contains classes and procedures for computing various statistical quantities related to the Geometric distribution.
!>  pm_distGeomCyclic       | This module contains classes and procedures for computing various statistical quantities related to the Cyclic Geometric distribution.
!>  pm_distKolm             | This module contains classes and procedures for computing various statistical quantities related to the Kolmogorov distribution.
!>  pm_distLogNorm          | This module contains classes and procedures for computing various statistical quantities related to the Lognormal distribution.
!>  pm_distLogUnif          | This module contains classes and procedures for computing various statistical quantities related to the LogUniform (or Reciprocal) distribution.
!>  pm_distMultiNorm        | This module contains classes and procedures for computing various statistical quantities related to the MultiVariate Normal (MVN) distribution.
!>  pm_distNegExp           | This module contains classes and procedures for computing various statistical quantities related to the Negative Exponential distribution.
!>  pm_distNorm             | This module contains classes and procedures for computing various statistical quantities related to the univariate Normal distribution.
!>  pm_distNormShell        | This module contains procedures and generic interfaces for computing the Multivariate Normal Shell density function or mixtures of such densities with varying parameters.
!>  pm_distPareto           | This module contains classes and procedures for computing various statistical quantities related to the (Truncated) Pareto distribution.
!>  pm_distPois             | This module contains classes and procedures for computing various statistical quantities related to the Poisson distribution.
!>  pm_distPower            | This module contains classes and procedures for computing various statistical quantities related to the (Truncated) Power distribution.
!>  pm_distPoweto           | This module contains classes and procedures for computing various statistical quantities related to the (Truncated) Power/Pareto distribution (hence the name Poweto).
!>  pm_distUnif             | This module contains classes and procedures for computing various statistical quantities related to the univariate Uniform distribution.
!>  pm_distUnifEll          | This module contains classes and procedures for computing various statistical quantities related to the MultiVariate Uniform Ellipsoid (MVUE) distribution.
!>  pm_distUnifPar          | This module contains classes and procedures for setting up and computing the properties of the MultiVariate Uniform Parallelepiped (MVUP) Distribution.
!>  pm_distUnifSphere       | This module contains classes and procedures for computing various statistical quantities related to the Uniform Spherical distribution.
!>  pm_ellipsoid            | This module contains classes and procedures for setting up and computing the properties of the hyper-ellipsoids in arbitrary dimensions.
!>  pm_err                  | This module contains classes and procedures for reporting and handling errors.
!>  pm_except               | This module contains procedures and generic interfaces and generic interfaces for testing for exceptional cases at runtime.
!>  pm_fftnr                | This module contains procedures and generic interfaces for computing the Discrete Fourier Transform of a real or complex sequence using radix-2 Cooley–Tukey Fast-Fourier Transform.
!>  pm_fftpack              | This module contains procedures and generic interfaces for computing the Discrete Fourier Transform of a real or complex sequence using a mixed-radix decimation-in-frequency Fast-Fourier Transform.
!>  pm_io                   | This module contains classes and procedures for input/output (IO) or generic display operations on standard displays or internal/external files.
!>  pm_kind                 | This module defines the relevant Fortran kind type-parameters frequently used in the ParaMonte library for the two standard supported Fortran and C-Fortran Interoperation (CFI) modes.
!>  pm_knn                  | This module contains procedures and generic interfaces for computing the nearest neighbor statistics of random samples.
!>  pm_lapack               | This module contains a set of generic interfaces to the LAPACK routines.
!>  pm_logicalCompare       | This module contains procedures and generic interfaces for performing a variety of logical comparison operations using logical values as if `.true.` evaluates to `1` and `.false.` evaluates to `0`.
!>  pm_math1mexp            | This module contains procedures and generic interfaces for computing `1 - exp(x)` more precisely for tiny `x`.
!>  pm_mathBeta             | This module contains classes and procedures for computing the mathematical Beta Function and its inverse.
!>  pm_mathCompare          | This module contains the procedures and interfaces for evaluating the relative or absolute proximity of two numeric values.
!>  pm_mathConst            | This module contains relevant mathematical constants.
!>  pm_mathCumPropExp       | This module contains the procedures and interfaces for computing the cumulative sum of the exponential of an array without undue numerical overflow.
!>  pm_mathCumSum           | This module contains the procedures and interfaces for computing the cumulative sum of an array
!>  pm_mathDivMul           | This module contains procedures and generic interfaces for evaluating the mathematical division and multiplication operators acting on `integer`, `complex`, or `real` values.
!>  pm_mathErf              | This module contains classes and procedures for computing the mathematical Inverse Error Function.
!>  pm_mathExp              | This module contains procedures and generic interfaces for computing the previous/next integer exponent for the given base that yields a number smaller/larger than the absolute input value.
!>  pm_mathFactorial        | This module contains procedures and generic interfaces for the Factorial function.
!>  pm_mathFactoring        | This module contains procedures and generic interfaces and generic interfaces for computing the prime factors of integers.
!>  pm_mathGamma            | This module contains procedures and generic interfaces for the Lower and Upper Incomplete Gamma functions.
!>  pm_mathLog1p            | This module contains procedures and generic interfaces for computing log(1 + x) more precisely for tiny `x`.
!>  pm_mathLogAddExp        | This module contains procedures and generic interfaces for adding two real or complex values without causing overflow or underflow.
!>  pm_mathLogSubExp        | This module contains procedures and generic interfaces for subtracting two real or complex values without causing overflow or underflow.
!>  pm_mathLogSumExp        | This module contains the procedures and interfaces for computing the natural logarithm of the sum of exponentials the elements of an array.
!>  pm_mathMinMax           | This module contains procedures and generic interfaces for finding the minimum and maximum of two input scalar values through lexical comparison.
!>  pm_mathNumSys           | This module contains procedures and generic interfaces for converting numbers to different bases in different numeral systems.
!>  pm_mathRoot             | This module contains classes and procedures for computing the roots of one-dimensional continuous mathematical functions using various root-finding methods.
!>  pm_mathRootTest         | This module contains a collection of example functions for testing or examining the root-finding routines of the ParaMonte library.
!>  pm_mathSqrt             | This module contains procedures and generic interfaces and generic interfaces for computing the square root of integers.
!>  pm_mathSubAdd           | This module contains procedures and generic interfaces for evaluating the mathematical operator ∓ acting on integer, complex, or real values.
!>  pm_mathUnsigned         | This module contains procedures and generic interfaces and generic interfaces for various operations with positive integers with results that have the same binary representation as an unsigned integer.
!>  pm_matrixChol           | This module contains procedures and generic interfaces for computing the Cholesky factorization of positive definite matrices.
!>  pm_matrixClass          | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.
!>  pm_matrixCopy           | This module contains procedures and generic interfaces relevant to copying (diagonal or upper/lower triangular) subsets of matrices of arbitrary intrinsic types and kinds from one matrix of arbitrary shape and packing format to another matrix of arbitrary shape and packing format.
!>  pm_matrixDet            | This module contains procedures and generic interfaces relevant to the computation of the determinants of square matrices.
!>  pm_matrixIndex          | This module contains procedures and generic interfaces for converting the indices of matrix elements between different packing and storage formats.
!>  pm_matrixInit           | This module contains procedures and generic interfaces relevant to generating and initializing matrices of arbitrary shapes `(:, :)`.
!>  pm_matrixInv            | This module contains abstract and concrete derived types and procedures related to the inversion of square matrices.
!>  pm_matrixLUP            | This module contains procedures and generic interfaces relevant to the partially LU Pivoted decomposition of matrix operations and linear algebra.
!>  pm_matrixMulAdd         | This module contains procedures and generic interfaces relevant to combined matrix-matrix or matrix-vector multiplication and addition.
!>  pm_matrixMulTri         | This module contains the procedures for multiplication of a square triangular matrix in various transpositions with a general matrix.
!>  pm_matrixPack           | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.
!>  pm_matrixSubset         | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.
!>  pm_matrixTrace          | This module contains procedures and generic interfaces for computing the additive or multiplicative trace of a given square matrix in arbitrary packing formats.
!>  pm_matrixTrans          | This module contains abstract and concrete derived types and procedures related to various common matrix transposition operations for which there is a corresponding matrix class defined in pm_matrixClass.
!>  pm_matrixUpdate         | This module contains procedures and generic interfaces relevant to arbitrary-rank updates to vectors, general matrices, or Symmetric/Hermitian triangular matrices of type integer, complex, and real of arbitrary type-kind parameters.
!>  pm_memory               | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for various search operations.
!>  pm_optimization         | This module contains procedures, generic interfaces, and types for numerical optimizations of mathematical functions.
!>  pm_option               | This module contains convenience functions for generating default values for optional arguments.
!>  pm_os                   | This module contains procedures and generic interfaces for inferring the processor operating system.
!>  pm_parallelism          | This module contains procedures and generic interfaces for facilitating parallel computations or computing the performance of the parallel Coarray/MPI/OpenMP algorithms.
!>  pm_paramonte            | This module contains procedures and data that provide general information about the ParaMonte library, its interfaces, and its build.
!>  pm_physUnit             | This module contains relevant physical constants.
!>  pm_polation             | This module contains procedures and data types for interpolation of finite samples of data.
!>  pm_polynomial           | This module contains procedures and generic interfaces for performing various mathematical operations involving polynomials.
!>  pm_quadPack             | This module contains classes and procedures for non-adaptive and adaptive global numerical quadrature and Cauchy Principal Value of 1D functions with various types of singularities and points of difficulties via the Gauss-Kronrod and Clenshaw-Curtis quadrature formulae.
!>  pm_quadRomb             | This module contains classes and procedures to perform numerical integrations.
!>  pm_quadTest             | This module contains a collection of interesting or challenging integrands for testing or examining the integration routines of the ParaMonte library.
!>  pm_sampleACT            | This module contains classes and procedures for computing properties related to the auto correlation time (ACT) of random sequences.
!>  pm_sampleAffinity       | This module contains classes and procedures for affine transformation of multivariate samples.
!>  pm_sampleCCF            | This module contains classes and procedures for computing properties related to the cross correlation of random samples.
!>  pm_sampleCor            | This module contains classes and procedures for computing properties related to the correlation matrices of random samples.
!>  pm_sampleCov            | This module contains classes and procedures for computing the properties related to the covariance matrices of a random sample.
!>  pm_sampleECDF           | This module contains classes and procedures for computing the Empirical Cumulative Distribution Function (ECDF) of an observational sample and the associated the various properties.
!>  pm_sampleMean           | This module contains classes and procedures for computing the first moment (i.e., the statistical mean) of random weighted samples.
!>  pm_sampleNorm           | This module contains classes and procedures for normalizing univariate or multivariate samples by arbitrary amounts along specific directions.
!>  pm_sampleQuan           | This module contains procedures and data types for computing sample quantile.
!>  pm_sampleScale          | This module contains classes and procedures for scaling (i.e., multiplying) univariate or multivariate samples by arbitrary amounts along specific directions.
!>  pm_sampleShift          | This module contains classes and procedures for shifting univariate or multivariate samples by arbitrary amounts along specific directions.
!>  pm_sampleVar            | This module contains classes and procedures for computing the properties related to the covariance matrices of a random sample.
!>  pm_sampleWeight         | This module contains the types, classes, and procedures relevant to weights of random samples.
!>  pm_sampling             | This module contains procedures and generic interfaces for the ParaMonte library sampler routines.
!>  pm_search               | This module contains abstract and concrete derived types that are required for compile-time resolution of procedures within the generic interfaces of the ParaMonte library for various search operations.
!>  pm_statest              | This module contains classes and procedures for performing various statistical tests.
!>  pm_str                  | This module contains classes and procedures for various string manipulations and inquiries.
!>  pm_strANSI              | This module contains procedures and generic interfaces for styling strings according for display on DEC VT100 or compatible terminals.
!>  pm_strASCII             | This module contains the uncommon and hardly representable ASCII characters as well as procedures for operating on strings that exclusively contain the 128 ASCII characters.
!>  pm_swap                 | This module contains procedures and generic interfaces for swapping values of intrinsic Fortran types of arbitrary kinds.
!>  pm_sysInfo              | This module contains procedures and generic interfaces for inferring the operating system kernel type, name, and other information.
!>  pm_sysPath              | This module contains classes and procedures for manipulating system file/folder paths.
!>  pm_sysShell             | This module contains procedures and generic interfaces for inferring the runtime system shell type and fetching information from the shell.
!>  pm_test                 | This module contains a simple unit-testing framework for the Fortran libraries, including the ParaMonte library.
!>  pm_timer                | This module contains the timer procedures and derived types to facilitate timing applications at runtime.
!>  pm_val2complex          | This module contains procedures and types for facilitating the conversion of values of different types (e.g., intrinsic Fortran string and logical) to complex values of different kinds.
!>  pm_val2int              | This module contains procedures and types for facilitating the conversion of values of different types (e.g., intrinsic Fortran string and logical) to integer values of different kinds.
!>  pm_val2logical          | This module contains procedures and types for facilitating the conversion of values of different types (e.g., intrinsic Fortran strings) to logical values of different kinds.
!>  pm_val2real             | This module contains procedures and types for facilitating the conversion of values of different types (e.g., intrinsic Fortran string and logical) to real values of different kinds.
!>  pm_val2str              | This module contains the generic procedures for converting values of different types and kinds to Fortran strings.
!>  pm_ziggurat             | This module contains procedures and generic interfaces for computing the Ziggurat set for for pseudo-random number sampling.
!>
!>  [⛓](#ParaMonteLangNamingConventions)
!>  \section ParaMonteLangNamingConventions ParaMonte Naming Conventions
!>
!>  +   The **CamelCase** naming style is enforced throughout the ParaMonte Fortran library.
!>
!>  [⛓⛓](#ParaMonteLangNamingConventionsVariables)
!>  \subsection ParaMonteLangNamingConventionsVariables ParaMonte Naming Conventions: Variables
!>
!>  The Fortran language is case-insensitive. However, by convention in this library,
!>
!>  +   All variables and procedure names except compile-time constants begin with an lowercase letter.
!>
!>  +   All constants and parameters should generally be typed in uppercase.
!>      <br>
!>      \example{constants}
!>      \code{.F90}
!>
!>          use iso_fortran_env, only: real64
!>          integer , parameter :: RK = real64
!>          real(RK), parameter :: PI = acos(-1._RK)
!>
!>      \endcode
!>
!>  +   The names of variables that **always** represent vectors of values **can**
!>      suffixed with `Vec` or `Vector` (for example, `proposalStd`, ...).<br>
!>
!>  +   The names of variables that **always** represent matrices of values **can**
!>      be suffixed with `mat` or `Matrix` (for example: `proposalCor`, ...).<br>
!>
!>  +   A significant attempt has been made to end all `logical` (Boolean) variables with a passive verb.<br>
!>      This is to ensure that the full variable name virtually forms a proposition.<br>
!>      In other words, a `logical` variable name should be an English-language statement that evaluates to either `.true.` or `.false.`.<br>
!>      For example, `parallelismMpiFinalizeEnabled` is one such proposition.<br>
!>      +   Occasionally, names that begin with the verb `is` can also be used to label `logical` objects.<br>
!>      +   But as a general rule, names that begin with a verb should be reserved for procedures.<br>
!>
!>  [⛓⛓](#ParaMonteLangNamingConventionsProcedures)
!>  \subsection ParaMonteLangNamingConventionsProcedures ParaMonte Naming Conventions: Procedures
!>
!>  +   Procedure (whether function, subroutine, or type-bound procedure) names
!>      should be descriptive of the action performed by the procedure. For example,
!>
!>      +   `getCov` means *generate (a) covariance matrix*.
!>      +   `setMatCopy` means *copy the matrix into the specified buffer*.
!>
!>  +   Procedure names should virtually always begin with a lowercase verb.
!>
!>  +   **Exceptions** to the above and below rules are allowed when,
!>      +   the procedure name is exceptionally famous, or
!>      +   it is very inconvenient to prefix the procedure name with a verb, or the prefixes `get`, or `set`.
!>
!>  [⛓⛓⛓](#ParaMonteLangNamingConventionsProceduresFunctions)
!>  \subsubsection ParaMonteLangNamingConventionsProceduresFunctions ParaMonte Naming Conventions: Functions
!>
!>  +   Function names should preferably begin with <b>`get`</b>.<br>
!>      The reasoning is simple: Functions *generate and obtain* a **new** object instead of changing (resetting) the state of an existing object.<br>
!>      For example, the function `getMatSym(mat) result(matsym)` generates a symmetric version of the input matrix and returns it as the function result.<br>
!>      **Exceptions** to this naming convention are allowed, for example, when a procedure is expected to exist only as a function (and not a subroutine) or when the function returns a `logical` result.<br>
!>
!>  +   Functions that return objects of type `logical` should be preferably prefixed with `is` or
!>      be named such that the name begins with a verb and reads as a proposition, evaluating to
!>      either `.true.` or `.false.`.
!>
!>  [⛓⛓⛓](#ParaMonteLangNamingConventionsProceduresSubroutines)
!>  \subsubsection ParaMonteLangNamingConventionsProceduresSubroutines ParaMonte Naming Conventions: Subroutines
!>
!>  +   The keyword <b>`get`</b> should be avoided as a prefix for subroutine names since unlike functions, subroutines do not
!>      generate and <i>get</i> a new object as their results, rather they <b>(re)set</b> the state of existing objects passed to them.<br>
!>  +   As such, <b>subroutine names should be always prefixed with `set`</b>, as in, for example, [setReplaced](@ref pm_arrayReplace::setReplaced).
!>
!>  [⛓](#ParaMonteLangAbbreviationGuidlines)
!>  \section ParaMonteLangAbbreviationGuidlines ParaMonte Abbreviation Guidelines
!>
!>  The following list of abbreviations is in alphabetical order to enable faster search:
!>
!>  +   The abbreviation `avg`      stands for **average** (rarely used).
!>  +   The abbreviation `cdf`      stands for **Cumulative Distribution Function** in the context of statistics. Example: `getNormCDF()`.
!>  +   The abbreviation `cho`      stands for **Cholesky factorization**. Example: `setChoLow()`.
!>  +   The abbreviation `chol`     stands for **Cholesky factorization**. Example: `setMatChol()`.
!>  +   The abbreviation `cor`      stands for **correlation**. Example: `getCor()`.
!>  +   The abbreviation `cov`      stands for **covariance**. Example: `getCov()`.
!>  +   The abbreviation `cum`      stands for **cumulative**. Example: `getCumSum()`.
!>  +   The abbreviation `coef`     stands for **coefficient**. Example: `corcoef_type()`.
!>  +   The abbreviation `def`      stands for **default** in variable names (mostly as a prefix `def_` or suffix `_def`).
!>  +   The abbreviation `def`      stands for **definite** (mostly in procedure names dealing with positive-definite matrices)
!>      \example{default}
!>      \code{.F90}
!>
!>          function getInvPosDefMat(mat) result(invPosDefMat)
!>              real, allocatable :: invPosDefMat(:,:) ! Def stands for definite.
!>          end function getInvPosDefMat
!>
!>          subroutine setAsserted(assertion, renabled)
!>              logical(LK), intent(in) :: assertion
!>              logical(LK), intent(in), optional :: renabled
!>              logical(LK) :: renabled_def ! def stands for default.
!>              renabled_def = .false._LK
!>              if (present(renabled)) renabled_def = renabled
!>          end subroutine
!>
!>      \endcode
!>  +   The abbreviation `den`      stands for **density**, mostly in the context of statistical procedures and objects. Example: `getLogProbDen()`.
!>  +   The abbreviation `det`      stands for **determinant**, mostly in the context of Matrix and linear algebra. Example: `getMatDet()`.
!>  +   The abbreviation `dia`      stands for **diagonal**, mostly in the context of matrix algebra, matrix packing, or Cholesky factorization. Example: `dia_type()`.
!>  +   The abbreviation `diag`     stands for **diagonal**, mostly as dummy argument in matrix algebra procedures.
!>  +   The abbreviation `desc`     stands for **description**, mostly as a dummy argument of [setAsserted](@ref pm_err::setAsserted) in tests.
!>  +   The abbreviation `diff`     stands for **difference**. Example: `setDisSortedExpDiff()`.
!>  +   The abbreviation `dist`     stands for **distance** or **distribution** depending on the context. Example: `DistMulti_type`.
!>  +   The abbreviation `eff`      stands for **effective**. Example: `effSamSize`.
!>  +   The abbreviation `exp`      stands for **exponential** or **exponentiated**. Example: `setDisSortedExpDiff()`.
!>  +   The abbreviation `hell`     stands for **Hellinger** in statistical distance computations. Example: `getDisHellSq()`.
!>  +   The abbreviation `herm`     stands for **hermitian** in matrix algebra.
!>  +   The abbreviation `ICE`      stands for **Internal Compiler Error**. It typically appears in the bug descriptions tagged via Doxygen command <tt>\\bug</tt>.
!>  +   The abbreviation `inv`      stands for **inverse**. Example: `getMatInv()`.
!>  +   The abbreviation `ks`       stands for **Kolmogorov-Smirnov** test. Example: `getProbKS()`.
!>  +   The abbreviation `lin`      stands for **linear**. Example: `getLinSpace()`.
!>  +   The abbreviation `low`      stands for **lower triangle of a matrix** or **lower limits**. Example: `setChoLow()`.
!>  +   The abbreviation `mahal`    stands for **Mahalanobis** distance. Example: `getDisMahalSq()`.
!>  +   The abbreviation `mat`      stands for **matrix**. Example: `getMatInv()`.
!>  +   The abbreviation `multi`    stands for **multivariate** mostly used in the context of statistical distributions. Example: `getMultiNormRand()`.
!>  +   The abbreviation `msn`      stands for **Multivariate Skew-Normal** mostly used in the context of the statistical MultiVariate Skew-Normal distribution.
!>  +   The abbreviation `mvn`      stands for **MultiVariate Normal** mostly used in the context of the statistical MultiVariate Normal distribution.
!>  +   The abbreviation `mvu`      stands for **MultiVariate Uniform** mostly used in the context of the statistical MultiVariate (ellipsoidal) Uniform distribution.
!>  +   The abbreviation `norm`     stands for **normal** in the context of statistical distributions or **normalization** factor. Example: `DistMultiNorm_type`.
!>  +   The abbreviation `normed`   stands for **normalized** mostly in the context of statistical samples. Example: `NormedSample`.
!>  +   The abbreviation `pdf`      stands for **Probability Density Function** in the context of statistics. Example: `getNormLogPDF()`.
!>  +   The abbreviation `pos`      stands for **positive**. Example: `getInvPosDefMat()`.
!>  +   The abbreviation `piwi`     stands for `piecewise`, mostly in the context of statistical applications. Example: `pm_PiwiPoweto()`.
!>  +   The abbreviation `prob`     stands for `probability`, mostly in the context of statistical applications. Example: `getLogProb()`.
!>  +   The abbreviation `proc`     stands for `procedure`, particularly, when it appears as the suffix `_proc` in `abstract interface` definitions.
!>      <br>
!>      \example{proc}
!>      \code{.F90}
!>
!>          use pm_paramonteInfo, only: getLogFunc_proc
!>          procedure(getLogFunc_proc) :: getLogFunc
!>
!>      \endcode
!>  +   The abbreviation `quan`     stands for **quantile**, mostly in the context of statistics. Example: `getParetoLogQuan()`.
!>  +   The abbreviation `rand`     stands for **random**, mostly in the context of statistics. Example: `getUnifRand()`.
!>  +   The abbreviation `ref`      stands for **reference**, mostly in the context of testings to represent the reference values for comparison. Example: `mean_ref`.
!>  +   The abbreviation `sam`      stands for **sample**, mostly in the context of statistics. Example: `effSamSize`.
!>  +   The abbreviation `sq`       stands for **squared**. Example: `getDisMahalSq()`.
!>  +   The abbreviation `stat`     stands for **statistics**. Example: `StatDRAM_type`.
!>  +   The abbreviation `std`      stands for **standard deviation**. Example: `StdVec`.
!>  +   The abbreviation `sym`      stands for **symmetric**.
!>  +   The abbreviation `symm`     stands for **symmetric**.
!>  +   The abbreviation `udf`      stands for **Unnormalized Density Function** in the context of statistics. Example: `getEggBoxLogUDF()`.
!>  +   The abbreviation `uni`      stands for **univariate**, mostly used in the context of statistical distributions. Example: `DistUni_type`.
!>  +   The abbreviation `unif`     stands for **uniform**, mostly in the context of the uniform statistical distribution. Example: `getUnifRand()`.
!>  +   The abbreviation `upp`      stands for **upper triangle of a matrix** or **upper limits**. Example: `setChoUpp()`.
!>  +   The abbreviation `vec`      stands for **vector**. Example: `stdVec`.
!>
!>  [⛓](#ParaMonteLangDeveloperWarnings)
!>  \section ParaMonteLangDeveloperWarnings ParaMonte Developer Guidelines and Warnings
!>
!>  The ParaMonte Fortran library development and guidelines are summarized in
!>  <a href="https://github.com/cdslaborg/paramonte/blob/main/src/fortran/CONTRIBUTING.md" target="_blank">CONTRIBUTING.md</a>.
!>
!>  [⛓](#ParaMonteLangDocumentationGuidelines)
!>  \section ParaMonteLangDocumentationGuidelines ParaMonte Fortran Documentation Guidelines
!>
!>  +   **Doxygen custom command orderings**.
!>
!>      +   The Doxygen tag `\brief` must always be the first line of the documentation of modules, types, and procedures.<br>
!>          Example: [pm_array](@ref pm_array).<br>
!>      +   The Doxygen tag `\details`, if it exists, must always immediately follow the Doxygen tag `\brief`.<br>
!>          Example: [pm_array](@ref pm_array).<br>
!>      +   The Doxygen tag `\param`, if any number of it exists, must always immediately follow the Doxygen tag `\brief` (or `\details` if it exists).<br>
!>          Example: [getMean()](@ref pm_sampleMean::getMean).<br>
!>      +   The Doxygen tag `\return`, must be exclusively used to indicate the return value of functions.<br>
!>          If it exists, it must appear immediately after the set of `\param` tags. Example: [getMean()](@ref pm_sampleMean::getMean).<br>
!>      +   If a generic interface is being documented, the ParaMonte custom command <tt>\\interface</tt> must appear immediately
!>          after the Doxygen `\return`, `\param`, `\details`, or `\brief` tags in the specified order, if any exists.<br>
!>      +   The Doxygen tag `\warning`, if any number of it exists, must immediately follow the Doxygen tag `\return` if it exists,
!>          otherwise `\param` if it exists, otherwise `\details` if it exists, otherwise `\brief`.<br>
!>          The `\warning` tag must be used to highlight situations that require special attention of the user,
!>          otherwise, there is a danger for the code section being documented to not behave normally as one may expect.<br>
!>      +   The Doxygen tag `\attention` has the same functionality and usage as `\warning`.<br>
!>          Therefore, `\warning` should be preferred wherever `\attention` is needed.<br>
!>          Exceptions are allowed and if they occur, the same documentation conventions as those of `\warning` also apply to the tag `\attention`.<br>
!>      +   The Doxygen tag `\remark`, if any number of it exists, must immediately follow the Doxygen tag `\warning` if it exists,
!>          otherwise the Doxygen tag `\return` if it exists, otherwise `\param` if it exists, otherwise `\details` if it exists, otherwise `\brief`.<br>
!>          The tag `\remark` should be reserved for explaining behavior that is **directly** related to the code segment being documented,
!>          but its knowledge is not so critical as warrant the use of a `\warning` tag.<br>
!>      +   The Doxygen tag `\note`, if it exists, must appear after all `\warning` and `\attention` and `\remark` tags and
!>          immediately before the ParaMonte custom command tag `\see` if it exists, otherwise immediately before <tt>\\example</tt> for examples (if it exists).<br>
!>      +   The Doxygen tag `\see`, if it exists, must appear after all `\warning` and `\remark` and `\note` tags.<br>
!>          If more than one item for the `\see` command exists, each must be written on a separate line and each line must end with the HTML line-break tag `<br>`.
!>          Example: See [below](#example-ParaMonteLangDocumentationGuidelines).<br>
!>      +   If any example exists, it must appear immediately after the `\see` tag, otherwise after `\note`, `\remark`, `\warning`, `\param`, `\details`, or `\brief` if any exists.<br>
!>          ParaMonte examples are initiated by the custom command <tt>\\example</tt> devised in the `config.txt` file of ParaMonte Doxygen documentation.<br>
!>          If the example exists in an external file, then it must be included via the Doxygen `\include` command, followed immediately by
!>          the ParaMonte custom Doxygen command <tt>\\compile</tt> which inserts the generic example compile commands for the example, followed optionally
!>          but immediately by the output file of the example inserted in the documentation via the `\include` command, followed immediately by the inclusion
!>          of any other visualization or postprocessing scripts and output.<br>
!>          <b>In all steps, it is imperative to not leave any empty lines between the successive commands of the example
!>          section</b>, designated by the <tt>\\example</tt>, otherwise, each empty line will start a new paragraph in the documentation.<br>
!>          Example: See [below](#example-ParaMonteLangDocumentationGuidelines).<br>
!>      +   The Doxygen `\test` tag, if any exists, must appear immediately after the example section designated by the <tt>\\example</tt> tag.<br>
!>      +   The Doxygen `\todo` tag, if any exists, must appear immediately after the `\test` tag or any other tag immediately preceding it.<br>
!>      +   The Doxygen `\bug` tag, if any exists, must appear immediately after the `\todo` tag or any other tag immediately preceding it.<br>
!>      +   The closing command of each documentation section must be the ParaMonte custom command <tt>\\final</tt> separated from the tags before and after by an empty line.<br>
!>      +   The Doxygen `\author` tag is the last command to appear in any documentation section, and it must preferably have the format exemplified in the example below.<br>
!>      <br>
!>      <br>
!>
!>  +   **ParaMonte Doxygen custom commands**.<br>
!>      To simplify documentation and avoid retyping certain frequently used keywords and sentences, a number of Doxygen aliases are predfined
!>      in the ParaMonte Doxygen `config.txt` file. These include (but are not limited to):
!>          +   <tt>\\warnpure</tt> Inserts a `\warning` about procedures that are `impure` when the library is built the preprocessor macro `CHECK_ENABLED=1`.
!>          +   <tt>\\elemental</tt> Inserts a `\remark` tag indicating that the procedure of interest is `elemental`.
!>          +   <tt>\\pure</tt> Inserts a `\remark` tag indicating that the procedure of interest is `pure`.
!>          +   <tt>\\interface</tt> Starts a **Possible calling interfaces** paragraph where different calling interfaces of a procedure can be listed.
!>          +   <tt>\\benchmark</tt> Starts a new **Benchmark** paragraph which is hyper-linked to the generic anchor `#benchmark` at the same location on the same page.
!>          +   <tt>\\benchmark{xxx}</tt> Starts a new **Benchmark** paragraph which is hyper-linked to the specific anchor `#benchmark-xxx` at the same location on the same page.
!>          +   <tt>\\benchmark{xxx, This is the benchmark title}</tt> Starts a new **Benchmark** paragraph which is hyper-linked to
!>              the specific anchor `#benchmark-xxx` at the same location on the same page with the title `This is the benchmark title`.
!>          +   <tt>\\example</tt> Starts a new **Example usage** paragraph which is hyper-linked to the generic anchor `#example` at the same location on the same page.
!>          +   <tt>\\example{xxx}</tt> Starts a new **Example usage** paragraph which is hyper-linked to the specific anchor `#example-xxx` at the same location on the same page.
!>          +   <tt>\\compile</tt> Inserts the set of example compile commands.
!>          +   <tt>\\output</tt> Inserts a title line for the output section of an example paragraph.
!>          +   <tt>\\postproc</tt> Inserts a title line for the postprocessing section of an example paragraph.
!>          +   <tt>\\abbr</tt> Inserts a `\remark` tag about the naming abbreviations used in the library.
!>          +   <tt>\\naming</tt> Inserts a `\remark` tag about the naming conventions used in the library.
!>          +   <tt>\\license</tt> Inserts a `\remark` tag about the generic licensing of the library.
!>          +   <tt>\\final</tt> Inserts the set of final generic remarks that should appear at the end of each documentation section.
!>          +   <tt>\\RK</tt>       Inserts a hyper-link reference \RK     to the default `real` kind used in the library.
!>          +   <tt>\\RK32</tt>     Inserts a hyper-link reference \RK32   to the `real32` real kind used in the library.
!>          +   <tt>\\RK64</tt>     Inserts a hyper-link reference \RK64   to the `real64` real kind used in the library.
!>          +   <tt>\\RK128</tt>    Inserts a hyper-link reference \RK128  to the `real128` real kind used in the library.
!>          +   <tt>\\CK</tt>       Inserts a hyper-link reference \CK     to the default `complex` kind used in the library.
!>          +   <tt>\\CK32</tt>     Inserts a hyper-link reference \CK32   to the `real32` complex kind used in the library.
!>          +   <tt>\\CK64</tt>     Inserts a hyper-link reference \CK64   to the `real64` complex kind used in the library.
!>          +   <tt>\\CK128</tt>    Inserts a hyper-link reference \CK128  to the `real128` complex kind used in the library.
!>          +   <tt>\\IK8</tt>      Inserts a hyper-link reference \IK8    to the `int8` integer kind used in the library.
!>          +   <tt>\\IK16</tt>     Inserts a hyper-link reference \IK16   to the `int16` integer kind used in the library.
!>          +   <tt>\\IK32</tt>     Inserts a hyper-link reference \IK32   to the `int32` integer kind used in the library.
!>          +   <tt>\\IK64</tt>     Inserts a hyper-link reference \IK64   to the `int64` integer kind used in the library.
!>          +   <tt>\\SKALL</tt>    Inserts a hyper-link reference to all major `character` kinds like: \SKALL.
!>          +   <tt>\\IKALL</tt>    Inserts a hyper-link reference to all major `integer`   kinds like: \IKALL.
!>          +   <tt>\\LKALL</tt>    Inserts a hyper-link reference to all major `logical`   kinds like: \LKALL.
!>          +   <tt>\\CKALL</tt>    Inserts a hyper-link reference to all major `complex`   kinds like: \CKALL.
!>          +   <tt>\\RKALL</tt>    Inserts a hyper-link reference to all major `real`      kinds like: \RKALL.
!>      <br>
!>      <br>
!>
!>      For an up-to-date list of all available aliases, check the value of the Doxygen `ALIASES` option in `config.txt` in the ParaMonte Fortran documentation repository.
!>
!>  +   **Escaping the Doxygen reserved characters**.<br>
!>      Doxygen has a set of reserved characters whose usage in the documentation must be handled properly.<br>
!>          +   Most importantly, the backslash character `\` begins a Doxygen command.<br>
!>              To print a backslash character to the output one should escape it via `\\`.<br>
!>          +   Also, the use of the percentage symbol `%` requires special care in some instances.<br>
!>              This is particularly important when defining Windows environment variables that should typically be enclosed with percentage character.<br>
!>              See the [pm_sysPath](@ref pm_sysPath) for instances of such definitions and how they are handled.<br>
!>      <br>
!>
!>      For more information, see [the relevant page on Doxygen documentation website](https://www.doxygen.nl/manual/commands.html#cmdfdollar).
!>      <br>
!>
!>  +   Avoid the insertion of an empty documentation line between any two lines of a single Doxygen paragraph.<br>
!>      This is crucial when the whole paragraph is indented by a vertical line as is done by Doxygen for `\warning`, `\remark`, `\note` and other similar tags.<br>
!>      \example{ParaMonteLangDocumentationGuidelines}
!>      The following is an example documentation for a procedure:
!>      \verbatim
!>
!>      !>  \brief
!>      !>  Generate and return the variance of the input array of shape `(np)` or `(nd,np)` or `(np,nd)` where `nd` is the number of
!>      !>  data dimensions (the number of data attributes) and `np` is the number of data points.
!>      !>
!>      !>  \param[in]  Sample  :   The input `contiguous` array of type `real` of kind \RKALL of shape `(np)`, `(nd,np)`, or `(np,nd)`
!>      !>                          containing the sample. If `Sample` is a 2D array, then the direction along which the variance is computed
!>      !>                          is dictated by the optional input argument `dim`.
!>      !>  \param[in]  Weight  :   The `contiguous` vector of shape `(np)` either type `real` of the same kind as the input `Sample` or type `integer`
!>      !>                          of kind \IKALL, containing the corresponding weight of each data points in `Sample`
!>      !>                          (**optional**, default = a vector of ones).
!>      !>  \param[in]  mean    :   The input scalar or `contiguous` vector of shape `(nd)` of the same type and kind as the input `Sample` containing
!>      !>                          the `Sample` mean along the (optionally) specified dimension `dim`. If the input `Sample` is a 1D array, then `mean`
!>      !>                          must be a scalar. Otherwise, if `mean` is a 2D array, then `mean` must be a vector whose size is the same as
!>      !>                          the size of at least one of the dimensions of `Sample`.
!>      !>                          (**optional**. If missing, then the input argument `shifted` must be present indicating whether the input `Sample`
!>      !>                          is already centered at the origin or it has to be shifted to the origin by the procedure).
!>      !>  \param[in]  shifted :   The input `logical` of default kind \LK indicating whether the input `Sample` is already centered at the origin or
!>      !>                          the it has to be shifted to the origin by the procedure (**optional**. If missing, then the input
!>      !>                          argument `mean` must be present).
!>      !>  \param[in]  biased  :   The input `logical` of default kind \LK indicating whether the output variance should be corrected for small sample-size
!>      !>                          bias. Set this argument to `.false.` to avoid biased variance computation, in particular, when the sample size `np`
!>      !>                          is small.
!>      !>  \param[in]  dim     :   An integer of default kind \IK indicating which dimension of the input `Sample` iterates over the individual data points.
!>      !>                          If `dim = 1` or `dim /= 2`, the input `Sample` is assumed to have the shape `(np,nd)`.
!>      !>                          If `dim = 2`, the input `Sample` is assumed to have the shape `(nd,np)`
!>      !>                          (**optional**, default = `2`. **This input argument is available only if the input `Sample` is a 2D array**.).
!>      !>
!>      !>  \return
!>      !>  `variance`          :   The output variance of the input sample of the same type and kind as the input `Sample`.
!>      !>                          It is a scalar only if the input `Sample` is a 1D array. Otherwise, it is an `allocatable` array of shape `(nd)`.
!>      !>
!>      !>  \warnpure
!>      !>
!>      !>  \note
!>      !>  One can also use the concise Fortran syntax to achieve the same goal as this function:
!>      !>  \code{.F90}
!>      !>
!>      !>       mean = sum(Weight*Sample) / sum(Weight)
!>      !>       variance = sum( (Weight*(Sample-mean))**2 ) / (sum(Weight)-1)
!>      !>
!>      !>  \endcode
!>      !>  But the above concise version will be slightly slower as it involves three loops instead of two.
!>      !>
!>      !>  \see
!>      !>  [getMean()](@ref pm_sampleMean::getMean)<br>
!>      !>
!>      !>  \example
!>      !>  \include{lineno} example/test_pm_sampleVar/getVar/main.F90
!>      !>  \compilef
!>      !>  \output
!>      !>  \include{lineno} example/test_pm_sampleVar/getVar/main.out.F90
!>      !>
!>      !>  \test
!>      !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
!>      !>
!>      !>  \todo
!>      !>  The performance of this code can improved.
!>      !>
!>      !>  \bug
!>      !>  This code used to have a well-known bug in version 1.1, but is now resolved.
!>      !>
!>      !>  \final
!>      !>
!>      !>  \author
!>      !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
!>
!>      \endverbatim
!>      <br>
!>      The above example documentation snippet will generate [an HTML similar to this documentation](@ref pm_sampleVar::getVar).<br>
!>      Note the lack of an empty line among the commands that immediately follow <tt>\\example</tt>.<br>
!>      This is essential to keep the entire example section in the same paragraph.
!>
!>  <br>
!>
!>  [⛓](#ParaMonteLangExamples)
!>  \section ParaMonteLangExamples ParaMonte Fortran Language Examples
!>
!>  The ParaMonte Fortran library ships with tens of thousands of example usage that are available in the `example/fortran` folder in the root directory of the project repository.<br>
!>  These examples are also available and discussed in the documentations of individual modules and procedures of this this documentation website.<br>
!>
!>  [⛓](#ParaMonteLangBenchmarks)
!>  \section ParaMonteLangBenchmarks ParaMonte Fortran Language Benchmarks
!>
!>  The ParaMonte Fortran library ships with a large number of performance benchmarks that are available in the `benchmark/fortran` folder in the root directory of the project repository.<br>
!>  These benchmarks are also available and discussed in the [benchmark listing page](@ref benchmarks) of this this documentation website.<br>
!>
!>  If you would like to see a relevant benchmark currently not included, [discuss it here](https://github.com/cdslaborg/paramonte/discussions)
!>  or [raise an issue here](https://github.com/cdslaborg/paramonte/issues) for consideration or volunteer to implement it!<br>
!>
!>  [⛓](#ParaMonteLangDocumentationTroubleshooting)
!>  \section ParaMonteLangDocumentationTroubleshooting ParaMonte Fortran Documentation Troubleshooting
!>
!>  <ol>
!>      <li>    **Side navigation pane disappears in some documentation pages.**<br>
!>              This issue most likely originates from the interference of browser addons with the documentation.<br>
!>              This issue is mostly observed on Firefox browsers.<br>
!>              If it occurs, open the page in *browser private mode* or use other (e.g., chrome-based) browsers.<br>
!>      <li>    **The ParaMonte Fortran documentation build fails because of a Doxygen lexer memory corruption leading to random segmentation faults.**<br>
!>              The current version of the ParaMonte library uses a [customized version](https://github.com/cdslaborg/doxygen) of Doxygen 1.9.3 documenter,
!>              specifically tailored to the needs of the Fortran interface of the ParaMonte library.<br>
!>              However, the Fortran lexer of this version of Doxygen is known to have a vicious memory corruption that leads to
!>              [random segmentation faults](https://github.com/doxygen/doxygen/issues/9298) when generating the documentation from the source files.<br>
!>              Although Doxygen developers have released newer versions of the documenter that is claimed to have resolved the bug,
!>              the bug appears to persist, at least within the ParaMonte library, even with the newer versions of the documenter.<br>
!>              **Resolution**: The random nature of the occurrence of such segfaults allows one to avoid the segfaults by slightly
!>              changing the documentation of source codes whenever it occurs.<br>
!>              For example, additions or removals as small as a line break \f$\ms{<br>}\f$ can resolve the bug.<br>
!>  </ol>
!>
!>  [⛓](#ParaMonteLangToDo)
!>  \section ParaMonteLangToDo ParaMonte Fortran ToDo List
!>
!>  For the full listing of all tasks to do see the dedicated [ToDo listing page](./todo.html).<br>
!>  The following are the library tasks that need to be accomplished.<br>
!>
!>  \todo
!>  \pvhigh
!>  The module `pm_distanceManhattan` for computing the Manhattan metric distance must be added to the library.<br>
!>  The module `pm_distanceMinkowski` for computing the Minkowski metric distance must be added to the library.<br>
!>
!>  \todo
!>  \pvhigh
!>  The module `pm_sampleConv` for timer series convolution must be added to the library.<br>
!>  The implementation of such module is straightforward and follows that of the existing module [pm_sampleCCF](@ref pm_sampleCCF).<br>
!>
!>  \todo
!>  \pvhigh
!>  The ParaNest and ParaDISE samplers must be added to the module [pm_sampling](@ref pm_sampling).<br>
!>  This is a task that only \AmirShahmoradi can complete.<br>
!>
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

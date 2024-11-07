!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains classes and procedures for computing various statistical quantities related to the univariate <b>Uniform distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the univariate <b>Uniform distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **continuous uniform distributions** or **rectangular distributions** are a family of symmetric probability distributions.<br>
!>  Such a distribution describes an experiment where there is an arbitrary outcome that lies between certain bounds.<br>
!>  The bounds are defined by the parameters, \f$a\f$ and \f$b\f$, which are the minimum and maximum values.<br>
!>  The interval can either be closed (i.e., \f$[a, b]\f$) or open (i.e., \f$(a, b)\f$).<br>
!>  Therefore, the distribution is often abbreviated as \f$U(a,b)\f$ where \f$U\f$ stands for uniform distribution.<br>
!>  The difference between the bounds defines the interval length.<br>
!>  All intervals of the same length on the distribution's support are equally probable.<br>
!>
!>  \note
!>  The Uniform distribution is the maximum entropy probability distribution for a random variable \f$x\f$ under no constraint other than that it is contained in the distribution's support.<br>
!>
!>  **Probability density function (PDF)**<br>
!>
!>  The **PDF** of the continuous uniform distribution is,
!>  \f{equation}{
!>      f(x) =
!>      \begin{cases}
!>          \frac{1}{b - a} &   \text{for} a\leq x \leq b ~, \\
!>          0               &   \text{for} x < a ~ \text{or} ~ x > b ~.
!>      \end{cases}
!>  \f}
!>
!>  The values of \f$f(x)\f$ at the two boundaries \f$a\f$ and \f$b\f$ are usually unimportant, because they do not alter the value of \f$\int_c^d f(x) dx\f$
!>  over any interval \f$[c,d]\f$ nor of \f$\int_a^b x f(x) dx\f$ nor of any higher moment.<br>
!>  Sometimes they are chosen to be zero, and sometimes chosen to be \f$\frac{1}{b-a}\f$.<br>
!>  The latter is appropriate in the context of estimation by the method of maximum likelihood.<br>
!>  In the context of Fourier analysis, one may take the value of \f$f(a)\f$ or \f$f(b)\f$ to be \f$\frac{1}{2(b - a)}\f$
!>  because then the inverse transform of many integral transforms of this uniform function will yield back the function itself,
!>  rather than a function which is equal *almost everywhere*, i.e., except on a set of points with zero measure.<br>
!>  Also, it is consistent with the sign function, which has no such ambiguity.<br>
!>  Any probability density function integrates to \f$1\f$.<br>
!>  Thus, the PDF of the continuous uniform distribution is graphically portrayed as a rectangle where \f$b − a\f$ is the base length and \f$\frac{1}{b-a}\f$ is the height.<br>
!>  As the base length increases, the height (the density at any particular value within the distribution boundaries) decreases.<br>
!>  In terms of mean \f$\mu\f$ and variance \f$\sigma^{2}\f$ the probability density function of the continuous uniform distribution is,
!>  \f{equation}{
!>      f(x) =
!>      \begin{cases}
!>          \frac{1}{2\sigma\sqrt{3}}   &   \text{for} -\sigma\sqrt{3} \leq x - \mu \leq \sigma\sqrt{3} ~, \\
!>          0                           &   \text{otherwise}                                            ~.
!>      \end{cases}
!>  \f}
!>
!>  **Cumulative distribution function (CDF)**<br>
!>
!>  The **CDF** of the continuous uniform distribution is,<br>
!>  \f{equation}{
!>      F(x) =
!>      \begin{cases}
!>           0                  &   \text{for}  x < a ~,\\
!>          \frac{x - a}{b - a} &   \text{for}  a\leq x\leq b ~, \\
!>          1                   &   \text{for}  x > b ~.
!>      \end{cases}
!>  \f}
!>
!>  In terms of mean \f$\mu\f$ and variance \f$\sigma^{2}\f$, the cumulative distribution function of the continuous uniform distribution is,
!>  \f{equation}{
!>      F(x) =
!>      \begin{cases}
!>          0                                                           &   \text{for}  x - \mu < -\sigma\sqrt{3} ~, \\
!>          \frac{1}{2}\left(\frac{x - \mu}{\sigma\sqrt{3}} + 1\right)  &   \text{for} -\sigma\sqrt{3} \leq x - \mu < \sigma\sqrt{3} ~, \\
!>          1                                                           &   \text{for} x - \mu \geq \sigma\sqrt{3}
!>      \end{cases}
!>  \f}
!>
!>  \note
!>  -#  The procedures under the generic interface [getUnifCDF](@ref pm_distUnif::getUnifCDF) of this module
!>      are `elemental` functions that accept `optional` arguments of arbitrary ranks.<br>
!>      As such, the procedures offer great flexibility in coding.<br>
!>      However, the elemental nature of the procedures impacts their runtime performance negatively.<br>
!>      See the benchmarks below for more information.<br>
!>  -#  The procedures under the generic interface [setUnifCDF](@ref pm_distUnif::setUnifCDF)
!>      are subroutines that accept a limited range of specific arguments ranks.<br>
!>      As such, they offer much better runtime performance compared to [getUnifCDF](@ref pm_distUnif::getUnifCDF)
!>      but have significantly less flexibility.<br>
!>  -#  **Which procedures should I use?**<br>
!>      The `elemental` procedures appear to incur no performance penalty with scalar arguments.
!>      However, there appears to exist a runtime performance penalty of ~2-3 times more than the rank-specific
!>      routines for array arguments, comparable to 10-20 CPU cycles.<br>
!>      These penalties are due to the looping that occur in the `elemental` procedures for array arguments.<br>
!>      However, this `elemental` performance penalty is likely **insignificant** in most practical cases,
!>      unless the `elemental` procedures are to be called on the order of tens of billions of times in a program,
!>      in which, case, the over all performance penalty, as of 2022, appears to be on the order of a few minutes or less.<br>
!>  -#  Note that the following benchmarks represent the worst case scenarios.<br>
!>      There may be situations where the compiler could inline the `elemental` procedures
!>      and remove the overhead of repeatedly calling the `elemental` function.
!>
!>  **Inverse Cumulative distribution function (ICDF) or Quantile Function**<br>
!>
!>  The Quantile function of continuous Uniform distribution is given by,
!>  \f{equation}{
!>      F^{-1}(p) = a + p(b - a) \quad \text{for} 0 < p < 1 ~.
!>  \f}
!>
!>  In terms of mean \f$\mu\f$ and variance \f$\sigma^{2}\f$, the Quantile function of the continuous uniform distribution is,
!>  \f{equation}{
!>      F^{-1}(p) = \sigma\sqrt{3} (2p - 1) + \mu \quad \text{for} 0 \leq p \leq 1 ~.
!>  \f}
!>
!>  **Random Number Generation (RNG)**<br>
!>
!>  This module contains two generic functional and subroutine interfaces
!>  <ol>
!>      <li>    [getUnifRand](@ref pm_distUnif::getUnifRand)
!>      <li>    [setUnifRand](@ref pm_distUnif::setUnifRand)
!>  </ol>
!>  for generating uniformly distributed random values of all intrinsic types and kinds supported by
!>  the Fortran standard and the processor (`character`, `integer`, `logical`, `complex`, `real`).<br>
!>  The functional interface is merely a wrapper around the generic subroutine interface.<br>
!>
!>  This module also contains four random number generator (RNG) algorithms that can be specified via the corresponding types,
!>  <ol>
!>      <li>    [rngf_type](@ref pm_distUnif::rngf_type) (the default intrinsic Fortran uniform RNG via `random_number()`)
!>      <li>    [splitmix64_type](@ref pm_distUnif::splitmix64_type)
!>      <li>    [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type)
!>      <li>    [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)
!>  </ol>
!>
!>  <b>Usage</b><br>
!>  <ol>
!>      <li>    [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) is the recommended RNG for all serial and parallel tasks.<br>
!>      <li>    [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) is the recommended RNG for tasks
!>              that mostly require `logical` random values, although it can be used for random value generation of any type and kind.<br>
!>      <li>    [splitmix64_type](@ref pm_distUnif::splitmix64_type) is the recommended RNG for initializing other RNGs
!>              or for simple serial tasks, although it can be used for random value generation of any type and kind.<br>
!>      <li>    The default Fortran RNG [rngf_type](@ref pm_distUnif::rngf_type), although flexible to use and fast,
!>              will not generate deterministic results across different compilers.<br>
!>  </ol>
!>
!>  \benchmarks
!>
!>  \benchmark{getUnifCDF_vs_setUnifCDF_D0, The runtime performance of scalar [getUnifCDF](@ref pm_distUnif::getUnifCDF) vs. [setUnifCDF](@ref pm_distUnif::setUnifCDF) without bounds}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0/main.F90
!>  \compilefb{getUnifCDF_vs_setUnifCDF_D0}
!>  \postprocb{getUnifCDF_vs_setUnifCDF_D0}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0/main.py
!>  \visb{getUnifCDF_vs_setUnifCDF_D0}
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0/benchmark.getUnifCDF_vs_setUnifCDF_D0.runtime.png width=1000
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0/benchmark.getUnifCDF_vs_setUnifCDF_D0.runtime.ratio.png width=1000
!>  \moralb{getUnifCDF_vs_setUnifCDF_D0}
!>      -#  The procedures under the generic interface [getUnifCDF](@ref pm_distUnif::getUnifCDF) are `elemental` functions with `optional` arguments.
!>          In the absence of the `optional` arguments, the default values are used, but the associated computations will be redundant.<br>
!>          As such, one expects the [getUnifCDF](@ref pm_distUnif::getUnifCDF) to perform less efficiently than
!>          the procedures under the generic interface [setUnifCDF](@ref pm_distUnif::setUnifCDF)
!>          which are rank-specific and carefully designed to avoid redundant computations.
!>
!>  \benchmark{getUnifCDF_vs_setUnifCDF_D1, The runtime performance of array [getUnifCDF](@ref pm_distUnif::getUnifCDF) vs. [setUnifCDF](@ref pm_distUnif::setUnifCDF) without bounds}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1/main.F90
!>  \compilefb{getUnifCDF_vs_setUnifCDF_D1}
!>  \postprocb{getUnifCDF_vs_setUnifCDF_D1}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1/main.py
!>  \visb{getUnifCDF_vs_setUnifCDF_D1}
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1/benchmark.getUnifCDF_vs_setUnifCDF_D1.runtime.png width=1000
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1/benchmark.getUnifCDF_vs_setUnifCDF_D1.runtime.ratio.png width=1000
!>  \moralb{getUnifCDF_vs_setUnifCDF_D1}
!>      -#  The procedures under the generic interface [getUnifCDF](@ref pm_distUnif::getUnifCDF) are `elemental` functions with `optional` arguments.<br>
!>          In the absence of the `optional` arguments, the default values are used, but the associated computations will be redundant.<br>
!>          Furthermore, `elemental` functions incur a performance penalty for input array arguments due to internal looping performed
!>          by the compiler to call the function repeatedly for different array elements.<br>
!>          As such, one expects the [getUnifCDF](@ref pm_distUnif::getUnifCDF) to perform less efficiently than
!>          the procedures under the generic interface [setUnifCDF](@ref pm_distUnif::setUnifCDF)
!>          which are rank-specific and carefully designed to avoid redundant computations.<br>
!>
!>  \benchmark{getUnifCDF_vs_setUnifCDF_D0_D0, The runtime performance of scalar [getUnifCDF](@ref pm_distUnif::getUnifCDF) vs. [setUnifCDF](@ref pm_distUnif::setUnifCDF) with bounds}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0_D0/main.F90
!>  \compilefb{getUnifCDF_vs_setUnifCDF_D0_D0}
!>  \postprocb{getUnifCDF_vs_setUnifCDF_D0_D0}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0_D0/main.py
!>  \visb{getUnifCDF_vs_setUnifCDF_D0_D0}
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0_D0/benchmark.getUnifCDF_vs_setUnifCDF_D0_D0.runtime.png width=1000
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D0_D0/benchmark.getUnifCDF_vs_setUnifCDF_D0_D0.runtime.ratio.png width=1000
!>  \moralb{getUnifCDF_vs_setUnifCDF_D0_D0}
!>      -#  The procedures under the generic interface [getUnifCDF](@ref pm_distUnif::getUnifCDF) are `elemental` functions with `optional` arguments.<br>
!>          In the presence of the `optional` arguments, the user-specified values are used.<br>
!>          Therefore, the costs of computations in [getUnifCDF](@ref pm_distUnif::getUnifCDF) are more comparable to
!>          the procedures under the generic interface [setUnifCDF](@ref pm_distUnif::setUnifCDF)
!>          which are rank-specific and carefully designed to avoid redundant computations.<br>
!>
!>  \benchmark{getUnifCDF_vs_setUnifCDF_D1_D0, The runtime performance of array [getUnifCDF](@ref pm_distUnif::getUnifCDF) vs. [setUnifCDF](@ref pm_distUnif::setUnifCDF) with bounds}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1_D0/main.F90
!>  \compilefb{getUnifCDF_vs_setUnifCDF_D1_D0}
!>  \postprocb{getUnifCDF_vs_setUnifCDF_D1_D0}
!>  \include{lineno} benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1_D0/main.py
!>  \visb{getUnifCDF_vs_setUnifCDF_D1_D0}
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1_D0/benchmark.getUnifCDF_vs_setUnifCDF_D1_D0.runtime.png width=1000
!>  \image html benchmark/pm_distUnif/getUnifCDF_vs_setUnifCDF_D1_D0/benchmark.getUnifCDF_vs_setUnifCDF_D1_D0.runtime.ratio.png width=1000
!>  \moralb{getUnifCDF_vs_setUnifCDF_D1_D0}
!>      -#  The procedures under the generic interface [getUnifCDF](@ref pm_distUnif::getUnifCDF) are `elemental` functions with `optional` arguments.<br>
!>          In the presence of the `optional` arguments, the user-specified values are used.<br>
!>          Although  the `elemental` functions incur a performance penalty for input array arguments due to internal looping performed
!>          by the compiler to call the function repeatedly for different array elements, the costs of computations in
!>          [getUnifCDF](@ref pm_distUnif::getUnifCDF) become more comparable to the procedures
!>          under the generic interface [setUnifCDF](@ref pm_distUnif::setUnifCDF) which
!>          are rank-specific and carefully designed to avoid redundant computations.<br>
!>
!>  \benchmark{setUnifRand_vs_random_number, The runtime performance of [setUnifRand](@ref pm_distUnif::setUnifRand) vs. Fortran intrinsic `random_number()`.}
!>  \include{lineno} benchmark/pm_distUnif/setUnifRand_vs_random_number/main.F90
!>  \compilefb{setUnifRand_vs_random_number}
!>  \postprocb{setUnifRand_vs_random_number}
!>  \include{lineno} benchmark/pm_distUnif/setUnifRand_vs_random_number/main.py
!>  \visb{setUnifRand_vs_random_number}
!>  \image html benchmark/pm_distUnif/setUnifRand_vs_random_number/benchmark.setUnifRand_vs_random_number.runtime.png width=1000
!>  \image html benchmark/pm_distUnif/setUnifRand_vs_random_number/benchmark.setUnifRand_vs_random_number.runtime.ratio.png width=1000
!>  \moralb{setUnifRand_vs_random_number}
!>      -#  The default RNG in the procedures under the generic interface [setUnifRand](@ref pm_distUnif::setUnifRand)
!>          are simply wrappers around the intrinsic random number generator of Fortran `random_number()`.<br>
!>          As such, [setUnifRand](@ref pm_distUnif::setUnifRand) for generating random real
!>          numbers has \f$\ms{5-10%}\f$ overhead with respect to the intrinsic `random_number()`.<br>
!>          Note Fortran does not have `integer`, `logical`, `complex`, or `character` uniform random RNG
!>          whereas [setUnifRand](@ref pm_distUnif::setUnifRand) provides a unified API for random numbers of all types.<br>
!>      -#  The RNGX is an acronym for [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) in the procedures under the generic interface [setUnifRand](@ref pm_distUnif::setUnifRand).<br>
!>          This random number generator, although unsafe for cryptographic purposes, is quite competitive and performant, even compared to the intrinsic Fortran compiler RNGs.<br>
!>
!>  \benchmark{splitmix64_type_vs_xoshiro256ss_type, The runtime performance of intrinsic `random_number()` vs. [splitmix64_type](@ref pm_distUnif::splitmix64_type) vs. [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).}
!>  \include{lineno} benchmark/pm_distUnif/splitmix64_type_vs_xoshiro256ss_type/main.F90
!>  \compilefb{splitmix64_type_vs_xoshiro256ss_type}
!>  \postprocb{splitmix64_type_vs_xoshiro256ss_type}
!>  \include{lineno} benchmark/pm_distUnif/splitmix64_type_vs_xoshiro256ss_type/main.py
!>  \visb{splitmix64_type_vs_xoshiro256ss_type}
!>  \image html benchmark/pm_distUnif/splitmix64_type_vs_xoshiro256ss_type/benchmark.splitmix64_type_vs_xoshiro256ss_type.runtime.png width=1000
!>  \image html benchmark/pm_distUnif/splitmix64_type_vs_xoshiro256ss_type/benchmark.splitmix64_type_vs_xoshiro256ss_type.runtime.ratio.png width=1000
!>  \moralb{splitmix64_type_vs_xoshiro256ss_type}
!>      -#  The [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) RNG
!>          greedily attempts to use as many randomly generated bits as possible in the output random values.<br>
!>      -#  The [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) RNG
!>          takes a wasteful approach of using at least one or more chunks of 64 randomly generated bits in the output random values.<br>
!>      -#  This fundamental difference between the two RNG types generally leads to faster random logical value generations with
!>          the greedy approach, because 64bits chunks translate to 64 logical values without updating the RNG state.<br>
!>      -#  However, the greedy approach leads to generally slower runtimes for real random value generation.<br>
!>      -#  Both greedy and wasteful RNGs appear to be much faster than the ParaMonte
!>          library wrappers for the implementations offered by \gfortran and \ifort.<br>
!>      -#  <b>Moral</b>: If your application requires many `logical` random number generation,
!>          use the greedy [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) RNG.<br>
!>          Conversely, if your application requires a mixture of random number generations of various types and kinds,
!>          use the wasteful [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) RNG.
!>
!>  \test
!>  [test_pm_distUnif](@ref test_pm_distUnif)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distUnif

    use pm_kind, only: SK, IK, LK, IK64

    implicit none

    !>  \cond excluded
    !private :: setStateNext, setStateJump
    !>  \endcond excluded

    character(*, SK), parameter :: MODULE_NAME = "@pm_distUnif"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Uniform
    !>  as defined in the description of [pm_distUnif](@ref pm_distUnif).
    !>
    !>  \details
    !>  See the documentation of [pm_distUnif](@ref pm_distUnif) for the definition of the Uniform distribution.
    !>
    !>  \interface{distUnif_type}
    !>  \code{.F90}
    !>
    !>      use pm_distUnif, only: distUnif_type
    !>      type(distUnif_type) :: distUnif
    !>
    !>      distUnif = distUnif_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distUnif_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distUnif_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of a univariate Standard Uniform distribution or a
    !>  Uniform distribution with the specified support via `lower` and `upper` input arguments at the specified input values.
    !>
    !>  \param[in]  x       :   The input scalar or array of the same shape as other input array arguments, of either <br>
    !>                          <ol>
    !>                              <li>    If `x` is `integer`, the discrete Uniform distribution CDF with support `[lower, upper]` will be returned.<br>
    !>                              <li>    If `x` is `integer`, the output argument `CDF` must be of type `real` of kind \RK.<br>
    !>                              <li>    If `x` is `real`, the continuous Uniform distribution CDF with support `[lower, upper]` will be returned.<br>
    !>                              <li>    If `x` is `complex`, the two real and imaginary components of `CDF` will correspond to two independent distributions.<br>
    !>                          </ol>
    !>  \param[in]  lower   :   The input scalar or array of the same shape as other input array arguments, of the same type and kind as `x`,
    !>                          representing the lower bound of the Uniform distribution.<br>
    !>                          (**optional**, default = `0`. It must be present **if and only if** the input argument `upper` is also present.)
    !>  \param[in]  upper   :   The input scalar or array of the same shape as other input array arguments, of the same type and kind as `x`,
    !>                          representing the upper bound of the Uniform distribution.<br>
    !>                          (**optional**, default = `1`. It must be present **if and only if** the input argument `lower` is also present.)
    !>
    !>  \return
    !>  `cdf`               :   The output scalar or array of the same shape as the input array arguments, of either <br>
    !>                          <ol>
    !>                              <li>    type `real` of default kind \RK (if the input value `x` is `integer` of kind \IKALL), or<br>
    !>                              <li>    type `complex` of the same kind as `x` (if the input value `x` is of type `complex` of kind \CKALL), <br>
    !>                              <li>    type `real` of the same kind as `x` (if the input value `x` is of type `real` of kind \RKALL), <br>
    !>                          </ol>
    !>                          containing the CDF of the specified discrete or continuous Uniform distribution.
    !>
    !>  \interface{getUnifCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distUnif, only: getUnifCDF
    !>
    !>      cdf = getUnifCDF(x)
    !>      cdf = getUnifCDF(x, lower, upper)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>
    !>  \example{getUnifCDF}
    !>  \include{lineno} example/pm_distUnif/getUnifCDF/main.F90
    !>  \compilef{getUnifCDF}
    !>  \output{getUnifCDF}
    !>  \include{lineno} example/pm_distUnif/getUnifCDF/main.out.F90
    !>  \postproc{getUnifCDF}
    !>  \include{lineno} example/pm_distUnif/getUnifCDF/main.py
    !>  \vis{getUnifCDF}
    !>  \image html pm_distUnif/getUnifCDF/getUnifCDF.IK.png width=700
    !>  \image html pm_distUnif/getUnifCDF/getUnifCDF.RK.png width=700
    !>  \image html pm_distUnif/getUnifCDF/getUnifCDF.CK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to input string arguments to make it compatible with [setUnifRand](@ref pm_distUnif::setUnifRand).
    !>
    !>  \final{getUnifCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! LU

    interface getUnifCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getUnifCDF_LU_IK5(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getUnifCDF_LU_IK4(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getUnifCDF_LU_IK3(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getUnifCDF_LU_IK2(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getUnifCDF_LU_IK1(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getUnifCDF_LU_CK5(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
        complex(CKG)                                :: cdf
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getUnifCDF_LU_CK4(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
        complex(CKG)                                :: cdf
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getUnifCDF_LU_CK3(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
        complex(CKG)                                :: cdf
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getUnifCDF_LU_CK2(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
        complex(CKG)                                :: cdf
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getUnifCDF_LU_CK1(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
        complex(CKG)                                :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getUnifCDF_LU_RK5(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getUnifCDF_LU_RK4(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getUnifCDF_LU_RK3(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getUnifCDF_LU_RK2(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getUnifCDF_LU_RK1(x, lower, upper) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_LU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
        real(RKG)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DD

    interface getUnifCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getUnifCDF_DD_IK5(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK
        integer(IKG), intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getUnifCDF_DD_IK4(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK
        integer(IKG), intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getUnifCDF_DD_IK3(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK
        integer(IKG), intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getUnifCDF_DD_IK2(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK
        integer(IKG), intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getUnifCDF_DD_IK1(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK
        integer(IKG), intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getUnifCDF_DD_CK5(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: cdf
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getUnifCDF_DD_CK4(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: cdf
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getUnifCDF_DD_CK3(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: cdf
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getUnifCDF_DD_CK2(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: cdf
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getUnifCDF_DD_CK1(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getUnifCDF_DD_RK5(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getUnifCDF_DD_RK4(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getUnifCDF_DD_RK3(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getUnifCDF_DD_RK2(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getUnifCDF_DD_RK1(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCDF_DD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of a univariate Standard Uniform distribution or a Uniform distribution
    !>  with the specified support via `lower` and `upper` input arguments at the specified input values.
    !>
    !>  \param[out] cdf     :   The output scalar or `contiguous` array of rank `1` of either <br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL (if the input value `x` is of type `complex`) or, <br>
    !>                              <li>    type `real` of kind \RKALL (if the input value `x` is of type `integer` or `real`), <br>
    !>                          </ol>
    !>                          containing the CDF of the specified discrete or continuous Uniform distribution.
    !>  \param[in]  x       :   The input scalar or `contiguous` array of the same shape as `cdf`, containing the values at which the CDF must be computed.<br>
    !>                          If `x` is of type `integer`, the CDF of the discrete Uniform distribution with support `[lower, upper]` will be returned.<br>
    !>                          If `x` is of type `integer`, the output argument `CDF` must be of type `real` of kind \RKALL.<br>
    !>                          If `x` is of type `real`, the output argument `CDF` must have the same type, kind, and rank as `x`,
    !>                          and will contain the CDF of the continuous Uniform distribution with support `[lower, upper)`.<br>
    !>                          If `x` is of type `complex`, the output argument `CDF` must have the same type, kind, and rank as `x`.<br>
    !>                          If `x` is of type `complex`, the two real and imaginary components of `CDF` will correspond to two independent distributions.<br>
    !>  \param[in]  lower   :   The input scalar of the same type and kind as `x`, representing the lower bound of the Uniform distribution.<br>
    !>                          (**optional**, default = `0`. If present, then `upper` must also be present.)
    !>  \param[in]  upper   :   The input scalar of the same type and kind as `x`, representing the upper bound of the Uniform distribution.<br>
    !>                          (**optional**, default = `1`. If present, then `lower` must also be present.)
    !>
    !>  \interface{setUnifCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distUnif, only: setUnifCDF
    !>
    !>      call setUnifCDF(cdf, x)
    !>      call setUnifCDF(cdf, x, lower, upper)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>
    !>  \example{setUnifCDF}
    !>  \include{lineno} example/pm_distUnif/setUnifCDF/main.F90
    !>  \compilef{setUnifCDF}
    !>  \output{setUnifCDF}
    !>  \include{lineno} example/pm_distUnif/setUnifCDF/main.out.F90
    !>  \postproc{setUnifCDF}
    !>  \include{lineno} example/pm_distUnif/setUnifCDF/main.py
    !>  \vis{setUnifCDF}
    !>  \image html pm_distUnif/setUnifCDF/setUnifCDF.IK.png width=700
    !>  \image html pm_distUnif/setUnifCDF/setUnifCDF.RK.png width=700
    !>  \image html pm_distUnif/setUnifCDF/setUnifCDF.CK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to input arguments with ranks higher than `1`.
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to input string arguments to make it compatible with [setUnifRand](@ref pm_distUnif::setUnifRand).
    !>
    !>  \final{setUnifCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! default range.

    interface setUnifCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK5_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK5_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK4_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK4_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK3_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK3_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK2_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK2_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK1_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK1_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK5_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK5_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK4_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK4_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK3_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK3_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK2_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK2_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK1_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK1_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK5_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK5_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK4_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK4_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK3_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK3_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK2_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK2_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK1_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK1_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK5_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK5_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK4_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK4_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK3_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK3_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK2_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK2_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK1_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK1_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK5_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK5_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK4_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK4_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK3_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK3_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK2_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK2_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK1_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK1_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_CK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_CK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_CK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_CK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_CK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D0_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK5_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK5_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK4_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK4_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK4_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK3_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK3_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK3_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK2_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK2_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK2_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK1_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK1_IK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK1_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK5_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK5_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK4_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK4_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK4_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK3_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK3_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK3_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK2_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK2_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK2_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK1_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK1_IK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK1_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK5_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK5_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK4_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK4_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK4_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK3_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK3_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK3_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK2_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK2_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK2_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK1_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK1_IK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK1_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK5_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK5_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK4_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK4_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK4_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK3_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK3_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK3_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK2_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK2_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK2_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK1_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK1_IK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK1_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK5_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK5_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK4_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK4_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK4_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK3_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK3_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK3_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK2_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK2_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK2_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK1_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK1_IK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK1_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_CK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_CK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_CK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_CK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_CK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifCDF_DD_D1_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_DD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! custom range.

    interface setUnifCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK5_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK5_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK4_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK4_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK3_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK3_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK2_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK2_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK1_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK1_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK5_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK5_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK4_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK4_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK3_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK3_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK2_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK2_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK1_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK1_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK5_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK5_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK4_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK4_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK3_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK3_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK2_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK2_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK1_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK1_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK5_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK5_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK4_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK4_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK3_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK3_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK2_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK2_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK1_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK1_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK5_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK5_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK4_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK4_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK3_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK3_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK2_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK2_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK1_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK1_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        integer(IKG), intent(in)                    :: x
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_CK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_CK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_CK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_CK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_CK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)                   :: cdf
        complex(CKG), intent(in)                    :: x
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D0_RK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK5_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK5_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK4_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK4_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK3_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK3_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK2_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK2_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK1_IK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK1_IK5
#endif
        use pm_kind, only: IKG => IK5, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK5_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK5_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK4_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK4_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK3_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK3_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK2_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK2_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK1_IK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK1_IK4
#endif
        use pm_kind, only: IKG => IK4, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK5_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK5_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK4_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK4_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK3_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK3_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK2_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK2_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK1_IK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK1_IK3
#endif
        use pm_kind, only: IKG => IK3, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK5_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK5_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK4_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK4_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK3_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK3_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK2_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK2_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK1_IK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK1_IK2
#endif
        use pm_kind, only: IKG => IK2, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK5_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK5_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK4_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK4_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK3_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK3_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK2_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK2_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED && IK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK1_IK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK1_IK1
#endif
        use pm_kind, only: IKG => IK1, RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        integer(IKG), intent(in)    , contiguous    :: x(:)
        integer(IKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_CK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_CK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_CK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_CK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_CK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)   , contiguous    :: cdf(:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG), intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK5(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK4(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK3(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK2(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifCDF_LU_D1_RK1(cdf, x, lower, upper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCDF_LU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: cdf(:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)   , intent(in)                    :: lower, upper
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `integer` of default kind containing
    !>  the number of binary digits of the `stream` component [Xoshiro256**](https://prng.di.unimi.it/) random number generator.<br>
    !>
    !>  \details
    !>  By definition, this number is `64`, because the type kind parameter of `stream` is \IK64.<br>
    !>
    !>  \see
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128)<br>
    !>  [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192)<br>
    !>  [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)<br>
    !>  [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)<br>
    !>
    !>  \final{xoshiro256ssStateSize}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    integer(IK)     , parameter :: xoshiro256ssStreamBitSize = int(bit_size(0_IK64), IK)

    !>  \brief
    !>  The constant scalar of type `integer` of default kind \IK containing
    !>  the size of the state vector of [Xoshiro256**](https://prng.di.unimi.it/) random number generator.<br>
    !>
    !>  \details
    !>  For more information see the documentation of
    !>  [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) and
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
    !>
    !>  \see
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128)<br>
    !>  [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192)<br>
    !>  [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)<br>
    !>  [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)<br>
    !>
    !>  \final{xoshiro256ssStateSize}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    integer(IK)     , parameter :: xoshiro256ssStateSize = 4_IK

    !>  \brief
    !>  The constant vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)
    !>  of type `integer` of kind \IK64 containing the state jump for the
    !>  [Xoshiro256**](https://prng.di.unimi.it/) random number generator.<br>
    !>
    !>  \details
    !>  This state jump can be passed to the constructor of [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) to request an RNG
    !>  whose state starts at `imageID * 2**128` steps (i.e., random number generations) ahead of the RNG constructed with `imageID = 1`.<br>
    !>  Using this jump, one can generate `2**128` independent RNG sequences each of which has a period of `2**128` in parallel applications.<br>
    !>  For more information see the documentation of
    !>  [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) and
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
    !>
    !>  The elements of this constant vector are obtained by transferring the following unsigned integers to signed values.<br>
    !>  \code{.F90}
    !>      integer(IK64)   , parameter :: xoshiro256ssJump128(xoshiro256ssStateSize) = [ transfer(Z"180ec6d33cfd0aba", 0_IK64) &
    !>                                                                                  , transfer(Z"d5a61266f0c9392c", 0_IK64) &
    !>                                                                                  , transfer(Z"a9582618e03fc9aa", 0_IK64) &
    !>                                                                                  , transfer(Z"39abdc4529b1661c", 0_IK64) ]
    !>  \endcode
    !>
    !>  \see
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128)<br>
    !>  [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192)<br>
    !>  [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)<br>
    !>  [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)<br>
    !>
    !>  \final{xoshiro256ssJump128}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    integer(IK64)   , parameter :: xoshiro256ssJump128(xoshiro256ssStateSize) = [ +1733541517147835066_IK64 &
                                                                                , -3051731464161248980_IK64 &
                                                                                , -6244198995065845334_IK64 &
                                                                                , +4155657270789760540_IK64 ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: xoshiro256ssJump128
#endif

    !>  \brief
    !>  The constant vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)
    !>  of type `integer` of kind \IK64 containing the state jump for the
    !>  [Xoshiro256**](https://prng.di.unimi.it/) random number generator.<br>
    !>
    !>  \details
    !>  This state jump can be passed to the constructor of [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) to request an RNG
    !>  whose state starts at `imageID * 2**192` steps (i.e., random number generations) ahead of the RNG constructed with `imageID = 1`.<br>
    !>  Using this jump, one can generate `2**64` independent RNG sequences each of which has a period of `2**192` in parallel applications.<br>
    !>  For more information see the documentation of
    !>  [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) and
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
    !>
    !>  The elements of this constant vector are obtained by transferring the following unsigned integers to signed values.<br>
    !>  \code{.F90}
    !>      integer(IK64)   , parameter :: xoshiro256ssJump192(xoshiro256ssStateSize) = [ transfer(Z"76e15d3efefdcbbf", 0_IK64) &
    !>                                                                                  , transfer(Z"c5004e441c522fb3", 0_IK64) &
    !>                                                                                  , transfer(Z"77710069854ee241", 0_IK64) &
    !>                                                                                  , transfer(Z"39109bb02acbe635", 0_IK64) ]
    !>  \endcode
    !>
    !>  \see
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128)<br>
    !>  [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192)<br>
    !>  [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)<br>
    !>  [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)<br>
    !>
    !>  \final{xoshiro256ssJump192}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    integer(IK64)   , parameter :: xoshiro256ssJump192(xoshiro256ssStateSize) = [ +8566230491382795199_IK64 &
                                                                                , -4251311993797857357_IK64 &
                                                                                , +8606660816089834049_IK64 &
                                                                                , +4111957640723818037_IK64 ]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` base derived type for defining various Uniform Random Number Generator (URNG) derived types.<br>
    !>
    !>  \see
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{rngu_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    type, abstract :: rngu_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances can be used to define/request
    !>  the default uniform random number generator (**RNG**) of the Fortran standard.<br>
    !>
    !>  \details
    !>  This type does not currently hold any components.<br>
    !>  It is merely used to signal/request the use of the intrinsic Fortran RNG.<br>
    !>  The Fortran programming language offers a single mechanism through the intrinsic procedure `random_seed(size = size, put = put, get = get)`
    !>  to set the seed of the intrinsic uniform random number generator of Fortran `random_number(harvest)`.<br>
    !>  However, the internal implementation of the random number generator can vary from compiler to another.<br>
    !>  Furthermore, the size and requirements for the initial state (seed) can be also different across compilers.<br>
    !>  As such, Fortran programming language also offers a new intrinsic procedure `random_init(repeatable, image_distinct)`.
    !>  <ol>
    !>      <li>    If `repeatable = .true.`, the seed is set to the same processor-dependent value
    !>              each time `random_init(repeatable = .true., ...)` is called from the same image.<br>
    !>              The program is, therefore, guaranteed to generate the same random number sequence after each call to `random_init()`.<br>
    !>              The term *same image* means a *single instance of program execution*.<br>
    !>              The sequence of random numbers is *different for repeated execution* of the program.<br>
    !>              If it is `.false.`, the seed is set to a processor-dependent value.<br>
    !>      <li>    If `image_distinct = .true.`, the seed is set to a processor-dependent value that is
    !>              distinct from th seed set by a call to `random_init()` in another (Coarray) image.<br>
    !>              If it is `image_distinct = .false.`, the seed is set to a value that does **not** depend which image called `random_init()`.<br>
    !>  </ol>
    !>  The Fortran intrinsic `random_init()` apparently does guarantee the generation
    !>  of the same random sequence in multiple independent executions of the program.<br>
    !>  This requires setting the RNG seed explicitly at the beginning of each program execution.<br>
    !>
    !>  The [rngf_type](@ref pm_distUnif::rngf_type) aims to facilitate the above
    !>  goal by offering a single interface that combines `random_init()` and `random_seed()` intrinsic routines.<br>
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of default kind \IK,
    !>                          containing a positive integer that serves as the starting point to generate the full RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a value determined by the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set identically on all images.)
    !>
    !>  \return
    !>   `rng`              :   The output scalar object of type [rngf_type](@ref pm_distUnif::rngf_type).<br>
    !>
    !>  \interface{rngf_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: rngf_type
    !>      type(rngf_type) :: rng
    !>
    !>      rng = rngf_type(seed = seed, imageID = imageID)
    !>      print *, getUnifRandState()
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This generic interface does not aim to replicate the behavior of the Fortran intrinsic `random_init()` but rather to extend it.<br>
    !>  If you wish to reset the RNG seed to a reproducible seed without specifying a scalar reference value for it via this generic interface,
    !>  simply call `random_init(repeatable = .true., ...)`.<br>
    !>  However, unlike `random_init()`, this generic interface can guarantee reproducibility of RNG sequences within
    !>  and between multiple program runs in serial or parallel Coarray/MPI/OpenMP, all compiled by the same compiler.<br>
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{rngf_type}
    !>  \include{lineno} example/pm_distUnif/rngf_type/main.F90
    !>  \compilef{rngf_type}
    !>  \output{rngf_type}
    !>  \include{lineno} example/pm_distUnif/rngf_type/main.out.F90
    !>
    !>  \final{rngf_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rngu_type) :: rngf_type
    end type

    !>  \brief
    !>  The scalar constant object of type [rngf_type](@ref pm_distUnif::rngf_type)
    !>  whose presence signified the use of the Fortran intrinsic random number generator (RNGF).<br>
    !>
    !>  \details
    !>  This constant is merely a convenience for making easier calls to routines that require a default RNGF.<br>
    !>
    !>  \interface{rngf}
    !>  \code{.F90}
    !>
    !>      use pm_distUnif, only: rngf, rngf_type
    !>      type(rngf_type) :: rng = rngf
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.8.0 20221119}
    !>  \desc
    !>  \ifort cannot handle the creation of a module constant of type [rngf_type](@ref pm_distUnif::rngf_type) as done for this object,
    !>  yielding the following error.<br>
    !>  \code{.sh}
    !>      error #9066: A generic function reference is not permitted in a constant expression.   [CONSTRUCTFRNG]
    !>  \endcode
    !>  GNU compiler compiles and runs the code without complaining.<br>
    !>  \remedy{2.0.0}
    !>  For now, the `parameter` attribute is removed from the declaration of [rngf](@ref pm_distUnif::rngf).<br>
    !>
    !>  \final{rngf}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(rngf_type) :: rngf! = rngf_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: rngf
#endif

    !>  \cond excluded
    interface rngf_type
        module procedure :: rngf_typer
    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar object of type [rngf_type](@ref pm_distUnif::rngf_type).
    !>
    !>  \details
    !>  This generic interface is the constructor of [rngf_type](@ref pm_distUnif::rngf_type).<br>
    !>  See the documentation of [rngf_type](@ref pm_distUnif::rngf_type) for example usage and interface.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{rngf_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface rngf_typer
    module function rngf_typer(seed, imageID) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: rngf_typer
#endif
        integer(IK) , intent(in), optional  :: seed, imageID
        type(rngf_type)                     :: self
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` base derived type for defining variants of
    !>  [Xoshiro256**](https://prng.di.unimi.it/) Uniform Random Number Generator derived types.<br>
    !>
    !>  \see
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{xoshiro256ss_type, The runtime performance of intrinsic `random_number()` vs. [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) vs. [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).}
    !>  \include{lineno} benchmark/pm_distUnif/xoshiro256ss_type/main.F90
    !>  \compilefb{xoshiro256ss_type}
    !>  \postprocb{xoshiro256ss_type}
    !>  \include{lineno} benchmark/pm_distUnif/xoshiro256ss_type/main.py
    !>  \visb{xoshiro256ss_type}
    !>  \image html benchmark/pm_distUnif/xoshiro256ss_type/benchmark.xoshiro256ss_type.runtime.png width=1000
    !>  \image html benchmark/pm_distUnif/xoshiro256ss_type/benchmark.xoshiro256ss_type.runtime.ratio.png width=1000
    !>  \moralb{xoshiro256ss_type}
    !>      -#  The [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) RNG
    !>          greedily attempts to use as many randomly generated bits as possible in the output random values.<br>
    !>      -#  The [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) RNG
    !>          takes a wasteful approach of using at least one or more chunks of 64 randomly generated bits in the output random values.<br>
    !>      -#  This fundamental difference between the two RNG types generally leads to faster random logical value generations with
    !>          the greedy approach, because 64bits chunks translate to 64 logical values without updating the RNG state.<br>
    !>      -#  However, the greedy approach leads to generally slower runtimes for real random value generation.<br>
    !>      -#  Both greedy and wasteful RNGs appear to be much faster than the ParaMonte
    !>          library wrappers for the implementations offered by \gfortran and \ifort.<br>
    !>      -#  <b>Moral</b>: If your application requires many `logical` random number generations,
    !>          use the greedy [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) RNG.
    !>          Conversely, if your application requires a mixture of random number generations of various types and kinds,
    !>          use the wasteful [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) RNG.
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{xoshiro256ss_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    type, abstract, extends(rngu_type) :: xoshiro256ss_type
        !>  \brief
        !>  The vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize)
        !>  of type `integer` of kind \IK64, containing the most recent RNG state.<br>
        integer(IK64)   :: state(xoshiro256ssStateSize) =   [ -5952639272145821898_IK64 &
                                                            , -2790978430781836137_IK64 &
                                                            , -4872796757339724681_IK64 &
                                                            , -6638731986642513151_IK64 ]
        !>  \brief
        !>  The scalar of type `integer` of kind \IK64, containing the most recently generated random 64-bit stream.
        integer(IK64)   :: stream
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for declaring and generating objects of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)
    !>  containing a unique instance of a [Xoshiro256**](https://prng.di.unimi.it/) random number generator (RNG).
    !>
    !>  \details
    !>  See also the documentation of [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer) for information on the constructor of this type.<br>
    !>
    !>  Xorshift random number generators, also called **shift-register generators**, are a class of pseudorandom
    !>  number generators that were invented by [George Marsaglia](https://en.wikipedia.org/wiki/George_Marsaglia).<br>
    !>  They are a subset of linear-feedback shift registers (LFSRs) which allow a particularly
    !>  efficient implementation in software without the excessive use of sparse polynomials.<br>
    !>  They generate the next number in their sequence by repeatedly taking the Exclusive OR (XOR)
    !>  of a number with a bit-shifted version of itself.<br>
    !>  This makes execution extremely efficient on modern computer architectures,
    !>  although it does not benefit efficiency in a hardware implementation.<br>
    !>  Like all LFSRs, the parameters have to be chosen very carefully in order to achieve a long period.<br>
    !>
    !>  The [xoshiro256**](https://prng.di.unimi.it/) RNG implemented in this derived type
    !>  is a subclass of Xorshift RNGs developed by David Blackman and Sebastiano Vigna.<br>
    !>  It is a 64-bit RNG that uses a carefully constructed linear transformation combined with shifts and rotations.<br>
    !>  This produces a computationally fast RNG with claimed excellent statistical qualities.<br>
    !>  Xoshiro256 has a period of \f$2^{256} - 1\f$ and and supports jumping the sequence in increments of \f$2^128\f$ and \f$2^192\f$.<br>
    !>  This allows the creation of many non-overlapping RNG subsequences for parallel applications.<br>
    !>
    !>  <b>xoshiro256** random seed</b>
    !>  The [xoshiro256**](https://prng.di.unimi.it/) state is determined by a vector of size
    !>  [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize) of type `integer` of kind \IK64.<br>
    !>  The RNG seed of [xoshiro256**](https://prng.di.unimi.it/) in this module is initialized based on either a user-specified scalar
    !>  initial seed value or an internal processor-dependent value based on the system clock.<br>
    !>  In either case, the seed is used as an input for another simple random number generator [splitMix64](@ref pm_distUnif::splitmix64_type)
    !>  and the output of this RNG is used as the initial state of [xoshiro256**](https://prng.di.unimi.it/).<br>
    !>
    !>  <b>Parallel applications</b>
    !>  [xoshiro256**](https://prng.di.unimi.it/) has a period of \f$2^{256} - 1\f$ and and supports jumping the sequence in increments of \f$2^128\f$ and \f$2^192\f$.<br>
    !>  This allows the creation of many non-overlapping RNG subsequences for parallel applications.<br>
    !>  [xoshiro256**](https://prng.di.unimi.it/) can be used in parallel applications by passing the `imageID` of the current processor to the RNG constructor.<br>
    !>  By default, the RNG constructor makes jumps of size `2**128` in the RRNG sequence to initialize subsequent RNGs on different processors.<br>
    !>  While such jump size is more than enough for any practical applications in modern days, the future galactic humans can request longer jumps
    !>  of size `2**192` by passing the non-default longer jump vector [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192) to the
    !>  constructor of the RNG.<br>
    !>  In either case, all RNGs on different processors must be initialized with the same original seed and jump vector (but with different processor
    !>  IDs set by the `imageID` argument to the RNG constructor) to ensure that the individual RNG sequences on different processors do not overlap.
    !>
    !>  <b>Usage instructions</b>
    !>  This derived type contains only the most recently updated state and random bit stream of the RNG.<br>
    !>  To generate random values of arbitrary intrinsic kinds (`character`, `integer`, `logical`, `complex`, `real`)
    !>  the user must,<br>
    !>  <ol>
    !>      <li>    declare an object of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) and
    !>              initialize the object via the type constructor [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)
    !>              (see below for the possible calling interfaces),
    !>      <li>    pass the generated RNG instance to the desired random number generating routines,
    !>              <ol>
    !>                  <li>    [getUnifRand](@ref pm_distUnif::getUnifRand),
    !>                  <li>    [setUnifRand](@ref pm_distUnif::setUnifRand),
    !>              </ol>
    !>              to generate the desired scalar or sequence of random values.<br>
    !>  </ol>
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of kind \IK64,
    !>                          containing an integer that serves as the starting point to generate the full deterministic RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible and deterministic,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a processor-dependent value based on the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set as if it is called in a serial application (or called on the first image).)
    !>  \param[in]  jump    :   The input vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize) of type `integer` of kind \IK64,
    !>                          whose value sets the jump size of the random number generator.<br>
    !>                          It can be,
    !>                          <ol>
    !>                              <li>    [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128), corresponding to a jump size of `imageID * 2**128`.<br>
    !>                                      This jump can be used to generate up to `2**128` unique RNG sequences in parallel, each with length `2**128`.<br>
    !>                              <li>    [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192), corresponding to a jump size of `imageID * 2**192`.<br>
    !>                                      This jump can be used to generate up to `2**64` unique RNG sequences in parallel, each with length `2**192`.<br>
    !>                          </ol>
    !>                          (**optional**. default = [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128))
    !>
    !>  \return
    !>   `rng`              :   The output scalar object (or array of objects, of the same rank and shape as the input array-like arguments)
    !>                          of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) containing an instance of a splitmix64
    !>                          random number generator.<br>
    !>
    !>  \interface{xoshiro256ssw_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: xoshiro256ssw_type
    !>      type(xoshiro256ssw_type) :: rng
    !>
    !>      rng = xoshiro256ssw_type(seed = seed, imageID = imageID, jump = jump)
    !>      rand = getUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  Although the components of this derived type are `public`, they are theoretically `protected`.<br>
    !>  The end users must not manipulate the component values at any stages of the random number generation.<br>
    !>
    !>  \note
    !>  Without initializing objects of this derived type, the generated RNGs will always be deterministic,
    !>  always yielding identical sequences.<br>
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{xoshiro256ssw_type}
    !>  \include{lineno} example/pm_distUnif/xoshiro256ssw_type/main.F90
    !>  \compilef{xoshiro256ssw_type}
    !>  \output{xoshiro256ssw_type}
    !>  \include{lineno} example/pm_distUnif/xoshiro256ssw_type/main.out.F90
    !>  \postproc{xoshiro256ssw_type}
    !>  \include{lineno} example/pm_distUnif/xoshiro256ssw_type/main.py
    !>  \vis{xoshiro256ssw_type}
    !>  \image html pm_distUnif/xoshiro256ssw_type/xoshiro256ssw_type.IK.png width=700
    !>  \image html pm_distUnif/xoshiro256ssw_type/xoshiro256ssw_type.CK.png width=700
    !>  \image html pm_distUnif/xoshiro256ssw_type/xoshiro256ssw_type.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \todo
    !>  \phigh
    !>  An illustration of the distribution of the probability of individual bits being `0` or `1` in the
    !>  mantissa of `real`-valued random numbers and `integer` random numbers must be added to the example.<br>
    !>
    !>  \final{xoshiro256ssw_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    type, extends(xoshiro256ss_type) :: xoshiro256ssw_type
    end type

    !>  \cond excluded
    interface xoshiro256ssw_type
        module procedure :: xoshiro256ssw_typer
    end interface
    !>  \endcond excluded

    !>  \brief
    !>  Generate, initialize, and return a scalar object of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).
    !>
    !>  \details
    !>  This generic interface is the constructor of [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
    !>  Upon return, the output object can be passed to [getUnifRand](@ref pm_distUnif::getUnifRand) and
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand) to generate uniformly-distributed
    !>  random values of various intrinsic types and kinds.<br>
    !>  See the documentation of [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) for example usage and calling interface.
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of kind \IK64,
    !>                          containing an integer that serves as the starting point to generate the full deterministic RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible and deterministic,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a processor-dependent value based on the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set as if it is called in a serial application (or called on the first image).)
    !>  \param[in]  jump    :   The input vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize) of type `integer` of kind \IK64,
    !>                          whose value sets the jump size of the random number generator.<br>
    !>                          It can be,
    !>                          <ol>
    !>                              <li>    [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128), corresponding to a jump size of `imageID * 2**128`.<br>
    !>                                      This jump can be used to generate up to `2**128` unique RNG sequences in parallel, each with length `2**128`.<br>
    !>                              <li>    [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192), corresponding to a jump size of `imageID * 2**192`.<br>
    !>                                      This jump can be used to generate up to `2**64` unique RNG sequences in parallel, each with length `2**192`.<br>
    !>                          </ol>
    !>                          (**optional**. default = [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128))
    !>
    !>  \return
    !>   `rng`              :   The output scalar object of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)
    !>                          containing an instance of a [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) random number generator.<br>
    !>
    !>  \interface{xoshiro256ssw_typer}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: xoshiro256ssw_type
    !>      use pm_distUnif, only: getUnifRand, setUnifRand
    !>      type(xoshiro256ssw_type) :: rng
    !>
    !>      rng = xoshiro256ssw_type(seed = seed, imageID = imageID, jump = jump)
    !>      rand = getUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>  [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)<br>
    !>  [rngf_typer](@ref pm_distUnif::rngf_typer)<br>
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)<br>
    !>
    !>  \final{xoshiro256ssw_typer}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    interface xoshiro256ssw_typer
    impure module function xoshiro256ssw_typer(seed, imageID, jump) result(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: xoshiro256ssw_typer
#endif
        use pm_kind, only: IKG => IK64
        integer(IKG)            , intent(in), optional  :: seed, jump(:)
        integer(IK)             , intent(in), optional  :: imageID
        type(xoshiro256ssw_type)                        :: rng
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for declaring and generating objects of type [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type)
    !>  containing a unique instance of a **greedy** [Xoshiro256**](https://prng.di.unimi.it/) random number generator (RNG).
    !>
    !>  \details
    !>  Unlike the `Xoshiro256**` algorithm as implemented by the derived type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>  the greedy version of the algorithm here does not waste any of the randomly generated 64 bits in each update the RNG state.<br>
    !>  See also the documentation of [xoshiro256ssg_typer](@ref pm_distUnif::xoshiro256ssg_typer) for information on the constructor of this type.<br>
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of kind \IK64,
    !>                          containing an integer that serves as the starting point to generate the full deterministic RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible and deterministic,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a processor-dependent value based on the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set as if it is called in a serial application (or called on the first image).)
    !>  \param[in]  jump    :   The input vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize) of type `integer` of kind \IK64,
    !>                          whose value sets the jump size of the random number generator.<br>
    !>                          It can be,
    !>                          <ol>
    !>                              <li>    [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128), corresponding to a jump size of `imageID * 2**128`.<br>
    !>                                      This jump can be used to generate up to `2**128` unique RNG sequences in parallel, each with length `2**128`.<br>
    !>                              <li>    [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192), corresponding to a jump size of `imageID * 2**192`.<br>
    !>                                      This jump can be used to generate up to `2**64` unique RNG sequences in parallel, each with length `2**192`.<br>
    !>                          </ol>
    !>                          (**optional**. default = [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128))
    !>
    !>  \return
    !>   `rng`              :   The output scalar object (or array of objects, of the same rank and shape as the input array-like arguments)
    !>                          of type [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) containing an instance of a splitmix64
    !>                          random number generator.<br>
    !>
    !>  \interface{xoshiro256ssg_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: xoshiro256ssg_type
    !>      type(xoshiro256ssg_type) :: rng
    !>
    !>      rng = xoshiro256ssg_type(seed = seed, imageID = imageID, jump = jump)
    !>      rand = getUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  Although the components of this derived type are `public`, they are theoretically `protected`.<br>
    !>  The end users must not manipulate the component values at any stages of the random number generation.<br>
    !>
    !>  \note
    !>  Without initializing objects of this derived type, the generated RNGs will always be deterministic,
    !>  always yielding identical sequences.<br>
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{xoshiro256ssg_type}
    !>  \include{lineno} example/pm_distUnif/xoshiro256ssg_type/main.F90
    !>  \compilef{xoshiro256ssg_type}
    !>  \output{xoshiro256ssg_type}
    !>  \include{lineno} example/pm_distUnif/xoshiro256ssg_type/main.out.F90
    !>  \postproc{xoshiro256ssg_type}
    !>  \include{lineno} example/pm_distUnif/xoshiro256ssg_type/main.py
    !>  \vis{xoshiro256ssg_type}
    !>  \image html pm_distUnif/xoshiro256ssg_type/xoshiro256ssg_type.IK.png width=700
    !>  \image html pm_distUnif/xoshiro256ssg_type/xoshiro256ssg_type.CK.png width=700
    !>  \image html pm_distUnif/xoshiro256ssg_type/xoshiro256ssg_type.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \todo
    !>  \phigh
    !>  An illustration of the distribution of the probability of individual bits being `0` or `1` in the
    !>  mantissa of `real`-valued random numbers and `integer` random numbers must be added to the example.<br>
    !>
    !>  \final{xoshiro256ssg_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    type, extends(xoshiro256ss_type) :: xoshiro256ssg_type
        !>  \brief
        !>  The scalar of type `integer` of default kind \IK, containing the
        !>  position of the first unused bit of the `stream` component of the RNG object.
        integer(IK)     :: pos = 0_IK
    end type

    !>  \cond excluded
    interface xoshiro256ssg_type
        module procedure :: xoshiro256ssg_typer
    end interface
    !>  \endcond excluded

    !>  \brief
    !>  Generate, initialize, and return a scalar object of type [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type).
    !>
    !>  \details
    !>  This generic interface is the constructor of [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type).<br>
    !>  Upon return, the output object can be passed to [getUnifRand](@ref pm_distUnif::getUnifRand) and
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand) to generate uniformly-distributed
    !>  random values of various intrinsic types and kinds.<br>
    !>  See the documentation of [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) for example usage and calling interface.
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of kind \IK64,
    !>                          containing an integer that serves as the starting point to generate the full deterministic RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible and deterministic,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a processor-dependent value based on the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set as if it is called in a serial application (or called on the first image).)
    !>  \param[in]  jump    :   The input vector of size [xoshiro256ssStateSize](@ref pm_distUnif::xoshiro256ssStateSize) of type `integer` of kind \IK64,
    !>                          whose value sets the jump size of the random number generator.<br>
    !>                          It can be,
    !>                          <ol>
    !>                              <li>    [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128), corresponding to a jump size of `imageID * 2**128`.<br>
    !>                                      This jump can be used to generate up to `2**128` unique RNG sequences in parallel, each with length `2**128`.<br>
    !>                              <li>    [xoshiro256ssJump192](@ref pm_distUnif::xoshiro256ssJump192), corresponding to a jump size of `imageID * 2**192`.<br>
    !>                                      This jump can be used to generate up to `2**64` unique RNG sequences in parallel, each with length `2**192`.<br>
    !>                          </ol>
    !>                          (**optional**. default = [xoshiro256ssJump128](@ref pm_distUnif::xoshiro256ssJump128))
    !>
    !>  \return
    !>   `rng`              :   The output scalar object of type [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type)
    !>                          containing an instance of a [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type) random number generator.<br>
    !>
    !>  \interface{xoshiro256ssg_typer}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: xoshiro256ssg_type
    !>      type(xoshiro256ssg_type) :: rng
    !>
    !>      rng = xoshiro256ssg_type(seed = seed, imageID = imageID, jump = jump)
    !>      rand = getUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand, lb, ub)
    !>      call setUnifRand(rng, rand)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>  [xoshiro256ssw_typer](@ref pm_distUnif::xoshiro256ssw_typer)<br>
    !>  [rngf_typer](@ref pm_distUnif::rngf_typer)<br>
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)<br>
    !>
    !>  \final{xoshiro256ssw_typer}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
    interface xoshiro256ssg_typer
    impure module function xoshiro256ssg_typer(seed, imageID, jump) result(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: xoshiro256ssg_typer
#endif
        use pm_kind, only: IKG => IK64
        integer(IKG)            , intent(in), optional  :: seed, jump(:)
        integer(IK)             , intent(in), optional  :: imageID
        type(xoshiro256ssg_type)                        :: rng
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for declaring and generating objects of type [splitmix64_type](@ref pm_distUnif::splitmix64_type)
    !>  containing a unique instance of an [splitmix64](https://doi.org/10.1145/2660193.2660195) random number generator (RNG).
    !>
    !>  \details
    !>  See also the documentation of [splitmix64_typer](@ref pm_distUnif::splitmix64_typer) for information on the constructor of this type.<br>
    !>
    !>  Splitmix64 is a pseudo-random number generator algorithm that originated from
    !>  the Java programming language and is used in many other programming languages.<br>
    !>  It is a fairly simple algorithm that is fast but unsuitable for cryptographic purposes and comparable tasks.<br>
    !>  The Splitmix64 algorithm is frequently used to calculate random initial states (seed) of other more complex pseudo-random number generators.<br>
    !>  The conventional splitmix64 holds one 64bit state (seed) variable and returns 64bits of random binary data upon subsequent calls.<br>
    !>  Splitmix64 is comparatively fast; It requires only 9 64-bit arithmetic/logical operations per 64 bits of random binary stream generation.<br>
    !>  A conventional linear RNG provides a generate method that returns one pseudorandom value and updates the state of the RNG.<br>
    !>  But the splitable RNG also has a second operation, split, that replaces the original RNG with two (seemingly) independent RNG.<br>
    !>  This is done by creating and returning a new such RNG and updating the state of the original object.<br>
    !>
    !>  <b>Applications</b>
    !>
    !>  Splitable RNG make it easy to organize the use of pseudorandom numbers in multithreaded programs structured using forkjoin parallelism.<br>
    !>  No locking or synchronization is required (other than the usual memory fence immediately after RNG creation).<br>
    !>  Because the generate method has no loops or conditionals, it is also suitable for SIMD or GPU implementation.<br>
    !>
    !>  <b>Usage instructions</b>
    !>  This derived type contains only the most recently updated state and random bit stream of the RNG.<br>
    !>  To generate random values of arbitrary intrinsic kinds (`character`, `integer`, `logical`, `complex`, `real`)
    !>  the user must,<br>
    !>  <ol>
    !>      <li>    declare an object of type [splitmix64_type](@ref pm_distUnif::splitmix64_type) and
    !>              initialize the object via the type constructor [splitmix64_typer](@ref pm_distUnif::splitmix64_typer)
    !>              (see below for the possible calling interfaces),
    !>      <li>    pass the generated RNG instance to the desired random number generating routines,
    !>              <ol>
    !>                  <li>    [getUnifRand](@ref pm_distUnif::getUnifRand),
    !>                  <li>    [setUnifRand](@ref pm_distUnif::setUnifRand),
    !>              </ol>
    !>              to generate the desired scalar or sequence of random values.<br>
    !>  </ol>
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of kind \IK64,
    !>                          containing an integer that serves as the starting point to generate the full RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible and deterministic,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a processor-dependent value based on the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set as if it is called in a serial application (or called on the first image).)
    !>
    !>  \return
    !>   `rng`              :   The output scalar object (or array of objects, of the same rank and shape as the input array-like arguments)
    !>                          of type [splitmix64_type](@ref pm_distUnif::splitmix64_type) containing an instance of a splitmix64
    !>                          random number generator.<br>
    !>
    !>  \interface{splitmix64_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: splitmix64_type
    !>      type(splitmix64_type) :: rng
    !>
    !>      rng = splitmix64_type(seed = seed, imageID = imageID)
    !>      print *, rng%state   ! The current RNG state.
    !>      print *, rng%stream  ! The current 64bit integer random number.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  Although the components of this derived type are `public`, they are theoretically `protected`.<br>
    !>  The end user must not change the values of the components of an object of type [splitmix64_type](@ref pm_distUnif::splitmix64_type)
    !>  at anytime during the random value generation.<br>
    !>
    !>  \warning
    !>  The Splitmix64 algorithm is not necessarily the best option for generating random values, particularly, for parallel applications.<br>
    !>  Use [xoshiro256**](@ref pm_distUnif::xoshiro256ssw_type) algorithm instead.<br>
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{splitmix64_type}
    !>  \include{lineno} example/pm_distUnif/splitmix64_type/main.F90
    !>  \compilef{splitmix64_type}
    !>  \output{splitmix64_type}
    !>  \include{lineno} example/pm_distUnif/splitmix64_type/main.out.F90
    !>  \postproc{splitmix64_type}
    !>  \include{lineno} example/pm_distUnif/splitmix64_type/main.py
    !>  \vis{splitmix64_type}
    !>  \image html pm_distUnif/splitmix64_type/splitmix64_type.IK.png width=700
    !>  \image html pm_distUnif/splitmix64_type/splitmix64_type.CK.png width=700
    !>  \image html pm_distUnif/splitmix64_type/splitmix64_type.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)<br>
    !>
    !>  \todo
    !>  \phigh
    !>  An illustration of the distribution of the probability of individual bits being `0` or `1` in the
    !>  mantissa of `real`-valued random numbers and `integer` random numbers must be added to the example.<br>
    !>
    !>  \final{splitmix64_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rngu_type) :: splitmix64_type
        !>  \brief
        !>  The scalar of type `integer` of kind \IK64, containing the most recently generated random 64-bit stream.
        integer(IK64)   :: stream
        !>  \brief
        !>  The scalar of type `integer` of kind \IK64, containing the most recent RNG state.<br>
        integer(IK64)   :: state = 324108011427370141_IK64
    end type

    !>  \cond excluded
    interface splitmix64_type
        module procedure :: splitmix64_typer
    end interface
    !>  \endcond excluded

    !>  \brief
    !>  Generate, initialize, and return a scalar object of type [splitmix64_type](@ref pm_distUnif::splitmix64_type).
    !>
    !>  \details
    !>  This generic interface is the constructor of [splitmix64_type](@ref pm_distUnif::splitmix64_type).<br>
    !>  Upon return, the output object can be passed to [getUnifRand](@ref pm_distUnif::getUnifRand) and
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand) to generate uniformly-distributed
    !>  random values of various intrinsic types and kinds using the [splitmix64](https://doi.org/10.1145/2660193.2660195) algorithm.<br>
    !>  See the documentation of [splitmix64_type](@ref pm_distUnif::splitmix64_type) for example usage and calling interface.
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of kind \IK64,
    !>                          containing an integer that serves as the starting point to generate the full RNG seed.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible and deterministic,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a processor-dependent value based on the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG seed uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG seed will be set as if it is called in a serial application (or called on the first image).)
    !>
    !>  \return
    !>   `rng`              :   The output scalar object (or array of objects, of the same rank and shape as the input array-like arguments)
    !>                          of type [splitmix64_type](@ref pm_distUnif::splitmix64_type) containing an instance of the splitmix64
    !>                          random number generator.<br>
    !>
    !>  \interface{splitmix64_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: splitmix64_type
    !>      type(splitmix64_type) :: rng
    !>
    !>      rng = splitmix64_type(seed = seed, imageID = imageID)
    !>      print *, rng%state   ! The current RNG state.
    !>      print *, rng%stream  ! The current 64bit integer random number.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{splitmix64_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface splitmix64_typer
    impure elemental module function splitmix64_typer(seed, imageID) result(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: splitmix64_typer
#endif
        use pm_kind, only: IKG => IK64
        integer(IKG)            , intent(in), optional  :: seed
        integer(IK)             , intent(in), optional  :: imageID
        type(splitmix64_type)                           :: rng
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    interface setStateNext

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setStateNextSM64(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateNextSM64
#endif
        use pm_kind, only: IK64
        type(splitmix64_type)   , intent(inout) :: rng
    end subroutine

    pure module subroutine setStateNextX256SSG(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateNextX256SSG
#endif
        use pm_kind, only: IK64
        type(xoshiro256ssg_type), intent(inout) :: rng
    end subroutine

    pure module subroutine setStateNextX256SSW(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateNextX256SSW
#endif
        use pm_kind, only: IK64
        type(xoshiro256ssw_type) , intent(inout) :: rng
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    interface setStateJump

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! default jump
    PURE module subroutine setStateJumpX256SSGDJ(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateJumpX256SSGDJ
#endif
        use pm_kind, only: IK64
        type(xoshiro256ssg_type), intent(inout) :: rng
    end subroutine

    ! arbitrary jump
    PURE module subroutine setStateJumpX256SSGAJ(rng, jump)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateJumpX256SSGAJ
#endif
        use pm_kind, only: IK64
        type(xoshiro256ssg_type), intent(inout) :: rng
        integer(IK64)           , intent(in)    :: jump(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! default jump
    PURE module subroutine setStateJumpX256SSWDJ(rng)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateJumpX256SSWDJ
#endif
        use pm_kind, only: IK64
        type(xoshiro256ssw_type), intent(inout) :: rng
    end subroutine

    ! arbitrary jump
    PURE module subroutine setStateJumpX256SSWAJ(rng, jump)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStateJumpX256SSWAJ
#endif
        use pm_kind, only: IK64
        type(xoshiro256ssw_type), intent(inout) :: rng
        integer(IK64)           , intent(in)    :: jump(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the size of the seed vector of the Fortran default random number generator (RNG).
    !>
    !>  \details
    !>  The functionality of this generic interface is equivalent to `random_seed(size = size)`.
    !>
    !>  \return
    !>   `seedSize` :   The output scalar of type `integer` of default kind \IK containing the size of the default Fortran RNG seed vector.
    !>
    !>  \interface{getUnifRandStateSize}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_distUnif, only: getUnifRandStateSize
    !>      integer(IK) :: seedSize
    !>
    !>      seedSize = getUnifRandStateSize()
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{getUnifRandStateSize}
    !>  \include{lineno} example/pm_distUnif/getUnifRandStateSize/main.F90
    !>  \compilef{getUnifRandStateSize}
    !>  \output{getUnifRandStateSize}
    !>  \include{lineno} example/pm_distUnif/getUnifRandStateSize/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{getUnifRandStateSize}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getUnifRandStateSize
    impure module function getUnifRandStateSizeDef() result(unifRandStateSize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandStateSizeDef
#endif
        use pm_kind, only: IK
        integer(IK)                 :: unifRandStateSize
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an `allocatable` array of rank `1` containing the state vector of the Fortran default random number generator (RNG) or,
    !>  optionally set the RNG state based on a reference input scalar seed, optionally distinctly on each processor.
    !>
    !>  \details
    !>  The procedures of this generic interface are merely wrappers around the procedures of the generic interface
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState) that additionally return the current RNG state.<br>
    !>  If both `seed` and `imageID` input arguments are missing, the procedures simply return the current RNG state.<br>
    !>  Otherwise, the behavior is identical to that of [setUnifRandState](@ref pm_distUnif::setUnifRandState).<br>
    !>
    !>  \param[in]  seed    :   See the documentation of the corresponding input argument to [setUnifRandState](@ref pm_distUnif::setUnifRandState).<br>
    !>  \param[in]  imageID :   See the documentation of the corresponding input argument to [setUnifRandState](@ref pm_distUnif::setUnifRandState).<br>
    !>
    !>  \return
    !>  `unifRandState`     :   The output `allocatable` vector of rank `1` whose length equals the `size` of the random state
    !>                          vector of the default Fortran RNG as returned by the Fortran intrinsic `random_seed(size = size)`.<br>
    !>
    !>  \interface{getUnifRandState}
    !>  \code{.F90}
    !>
    !>      use pm_lkind, only: IK
    !>      use pm_distUnif, only: getUnifRandState
    !>      integer(IK), allocatable :: unifRandState(:)
    !>
    !>      unifRandState = getUnifRandState(seed = seed, imageID = imageID)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{getUnifRandState}
    !>  \include{lineno} example/pm_distUnif/getUnifRandState/main.F90
    !>  \compilef{getUnifRandState}
    !>  \output{getUnifRandState}
    !>  \include{lineno} example/pm_distUnif/getUnifRandState/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{getUnifRandState}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getUnifRandState
    impure module function getUnifRandStateDef(seed, imageID) result(unifRandState)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandStateDef
#endif
        use pm_kind, only: IK, LK
        integer(IK) , intent(in)    , optional  :: seed, imageID
        integer(IK) , allocatable               :: unifRandState(:)
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the state of the Fortran default random number generator (RNG) to a random value or to an optionally deterministic,
    !>  optionally processor-dependent value based on the user-specified input scalar seed and processor ID.
    !>
    !>  \details
    !>  The procedures of this generic interface offer a convenient interface for the following tasks:
    !>  <ol>
    !>      <li>    When all optional arguments are missing, this generic interface
    !>              (re)sets the state of the default Fortran RNG to random values.<br>
    !>      <li>    When the optional argument `seed` is present and `imageID` is missing,
    !>              this generic interface sets the RNG state vector to a deterministic value
    !>              exclusively based on the input scalar `seed`.<br>
    !>              The generated state vector is set identically on all
    !>              images/processes/threads that individually call this generic interface.<br>
    !>      <li>    When both optional arguments `seed` and `imageID` are present,
    !>              this generic interface sets the RNG state vector to a deterministic value
    !>              unique to each process/thread/image (based on the specified input process-unique `imageID`).<br>
    !>  </ol>
    !>
    !>  \param[in]  seed    :   The input scalar of type `integer` of default kind \IK,
    !>                          containing a positive integer that serves as the starting point to generate the full RNG state.<br>
    !>                          Specify this input argument if you wish to make random simulations reproducible,
    !>                          even between multiple independent runs of the program compiled by the same compiler.<br>
    !>                          (**optional**. If missing, it is set to a value determined by the current date and time.)
    !>  \param[in]  imageID :   The input positive scalar `integer` of default kind \IK containing the ID of the current image/thread/process.<br>
    !>                          This can be,
    !>                          <ol>
    !>                              <li>    The Coarray image ID as returned by Fortran intrinsic `this_image()` within a global team of Coarray images.<br>
    !>                              <li>    The MPI rank of the processor (plus one) as returned by the MPI library intrinsic `mpi_comm_rank()`.<br>
    !>                              <li>    The OpenMP thread number as returned by the OpenMP library intrinsic `omp_get_thread_num()`.<br>
    !>                              <li>    Any (positive) integer that uniquely identifies the current processor from other processes.<br>
    !>                          </ol>
    !>                          The image/process/thread ID can be readily obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          This number will be used to set the RNG state uniquely on each processor.<br>
    !>                          (**optional**. If missing, the RNG state will be set identically on all images.)
    !>
    !>  \interface{setUnifRandState}
    !>  \code{.F90}
    !>
    !>      use pm_lkind, only: IK
    !>      use pm_distUnif, only: setUnifRandState
    !>
    !>      call setUnifRandState(seed = seed, imageID = imageID)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < imageID` must hold for the corresponding input arguments.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{setUnifRandState}
    !>  \include{lineno} example/pm_distUnif/setUnifRandState/main.F90
    !>  \compilef{setUnifRandState}
    !>  \output{setUnifRandState}
    !>  \include{lineno} example/pm_distUnif/setUnifRandState/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \final{setUnifRandState}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setUnifRandState
    impure module subroutine setUnifRandStateDef(seed, imageID)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandStateDef
#endif
        use pm_kind, only: IK
        integer(IK) , intent(in)    , optional  :: seed, imageID
    end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar or a `contiguous` array of rank `1` of length `s1` of randomly uniformly distributed
    !>  discrete `logical`, `integer`, `character` value(s), or continuous `real` or `complex value(s) within the specified input range.<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar of type,
    !>                              <ol>
    !>                                  <li>    [rngf_type](@ref pm_distUnif::rngf_type), or
    !>                                  <li>    [splitmix64_type](@ref pm_distUnif::splitmix64_type), or
    !>                                  <li>    [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type),
    !>                                  <li>    [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                              </ol>
    !>                              containing the user-specified random number generator algorithm to be used.<br>
    !>                              The user must initialize the object with the corresponding type constructors if non-deterministic RNG are desired.
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type))
    !>  \param[in]      lb      :   The input scalar (or array of the same shape as the desired output `rand` and other input array-like arguments)
    !>                              of the same type and kind as the output `rand`, representing the lower bound of the Uniform distribution.<br>
    !>                              If the output random value `rand` is to be of type `logical`, then `lb`, if present, must be `logical(.false., kind = kind(rand))`.<br>
    !>                              If the input argument `s1` is present, then `lb` must be a scalar.<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]      ub      :   The input scalar or array of the same type, kind, rank, and shape as `lb`, representing the upper bound of the Uniform distribution.<br>
    !>                              If the input argument `s1` is present, `ub` must be a scalar.<br>
    !>                              If the output random value `rand` is to be of type `logical`, then `ub`, if present, must be `logical(.true., kind = kind(rand))`.<br>
    !>                              (**optional**, default = `.true.`)
    !>  \param[in]      s1      :   The input scalar of type `integer` of default kind \IK, representing the size of the output `rand` array along its first dimension.<br>
    !>                              (**optional**. It must be present if `s2` is present. If missing, the rank and size of the output `rand` is that of `lb` or `ub` with a non-zero rank.)
    !>  \param[in]      s2      :   The input scalar of type `integer` of default kind \IK, representing the size of the output `rand` array along its second dimension.<br>
    !>                              (**optional**. It must be present if `s3` is present. If missing, the rank and size of the output `rand` is that of `lb` or `ub` with a non-zero rank.)
    !>
    !>  \return
    !>   `rand`                 :   The output scalar or vector of rank `1` of length `s1` or array of the same shape as the non-zero rank of the input `lb` and `ub`,
    !>                              containing the uniformly-distributed random output value(s).<br>
    !>                              If `lb` and `ub` are missing, then the output `rand` is of type `logical` of default kind \LK.<br>
    !>                              Otherwise, it is either of,
    !>                              <ol>
    !>                                  <li>    type `character` of the same kind as `lb` and `ub`, of the same length type parameter `len` as that of `lb` and `ub` or, <br>
    !>                                  <li>    type `integer` of the same kind as `lb` and `ub`, or<br>
    !>                                  <li>    type `logical` of the same kind as `lb` and `ub`, or<br>
    !>                                  <li>    type `complex` of the same kind as `lb` and `ub`, or<br>
    !>                                  <li>    type `real` of the same kind as `lb` and `ub`.<br>
    !>                              </ol>
    !>                              Note that,
    !>                              <ol>
    !>                                  <li>    If `lb` is of type `character`, then the random characters will be drawn from the **character collating sequence of the processor**.<br>
    !>                                  <li>    If `lb` and `ub` are of type `integer`, then `rand` will be in the interval `[lb, ub]`.<br>
    !>                                  <li>    If `lb` and `ub` are of type `real`, then `rand` will be in the interval `[lb, ub)`.<br>
    !>                                  <li>    If `lb` and `ub` are of type `character`, then `rand` will be in the interval `[lb, ub]` as defined by the processor collating sequence.<br>
    !>                                  <li>    If `lb` and `ub` are of type `complex`, then `rand%%re` and `rand%%im` will be in the intervals `[lb%%re, ub%%re)` and `[lb%%im, ub%%im)` respectively.<br>
    !>                                  <li>    If `lb` and `ub` are both missing, then the output will be a `logical` of default kind \LK with possible values `.false._LK` and `.true._LK`.
    !>                              </ol>
    !>
    !>  \interface{getUnifRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnif, only: getUnifRand
    !>
    !>      rand = getUnifRand() ! Non-elemental, random `logical` output scalar.
    !>      rand = getUnifRand(lb, ub) ! Elemental output of the same type and kind as `lb` and `ub`: character, complex, logical, integer, real.
    !>      rand(1:s1) = getUnifRand(lb, ub, s1) ! All input arguments must be scalars. Output is a vector of size `s1`, of the same type and kind as `lb` and `ub`.
    !>      rand(1:s1,1:s2) = getUnifRand(lb, ub, s1, s2) ! All input arguments must be scalars. Output is a matrix of size `(s1,s2)`, of the same type and kind as `lb` and `ub`.
    !>      rand(1:s1,1:s2,1:s3) = getUnifRand(lb, ub, s1, s2, s3) ! All input arguments must be scalars. Output is a cube of size `(s1,s2,s3)`, of the same type and kind as `lb` and `ub`.
    !>
    !>      rand = getUnifRand(rng) ! Non-elemental, random `logical` output scalar.
    !>      rand = getUnifRand(rng, lb, ub) ! Elemental output of the same type and kind as `lb` and `ub`: character, complex, logical, integer, real.
    !>      rand(1:s1) = getUnifRand(rng, lb, ub, s1) ! All input arguments must be scalars. Output is a vector of size `s1`, of the same type and kind as `lb` and `ub`.
    !>      rand(1:s1,1:s2) = getUnifRand(rng, lb, ub, s1, s2) ! All input arguments must be scalars. Output is a matrix of size `(s1,s2)`, of the same type and kind as `lb` and `ub`.
    !>      rand(1:s1,1:s2,1:s3) = getUnifRand(rng, lb, ub, s1, s2, s3) ! All input arguments must be scalars. Output is a cube of size `(s1,s2,s3)`, of the same type and kind as `lb` and `ub`.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `len(lb) == len(ub) .or. len(lb) == 1 .or. len(ub) == 1` for the corresponding input arguments of type `character`.<br>
    !>  The condition `lb <= ub` must hold for the corresponding input arguments where
    !>  logical values are compared by the procedures of module [pm_logicalCompare](@ref pm_logicalCompare) and
    !>  complex values are compared by the procedures of module [pm_complexCompareAll](@ref pm_complexCompareAll).<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>  The procedures of this generic interface are non-`elemental` when the argument `rng` is present.<br>
    !>
    !>  \remark
    !>  The procedures under this generic interface are carefully designed to avoid possible overflow due to
    !>  specifying huge negative and positive `lb` and `ub` limits of type `integer`, `complex`, `real`.<br>
    !>  This is possible at the cost of making the random number generation slightly more expensive (by a few CPU cycles, equivalent to and extra multiplication).<br>
    !>
    !>  \remark
    !>  The procedures under the generic interface are all `elemental` when the input arguments `lb` and `ub` are present and `s1` is missing.<br>
    !>
    !>  \remark
    !>  When the input `s1` argument is present, the procedures under this generic interface will be non-elemental,
    !>  and the input arguments `lb` and `ub` must be scalars, because the size of the output `rand` object is determined by `s1`.<br>
    !>  Additionally, when there is no input argument, the output will always be a scalar `logical` of default kind \LK.<br>
    !>
    !>  \note
    !>  The interface `getUnifRand()` corresponds to coin-flipping experiment with two possible outputs: `.true.` (head) and `.false.` (tail).<br>
    !>  See [isHead](@ref pm_distBern::isHead) for generating values from a **biased-coin** flipping experiment.<br>
    !>
    !>  \note
    !>  By default random characters are generated from the ASCII collating sequence to ensure portability across compilers and platforms.<br>
    !>  If random uniform characters from the processor collating sequence are desired, specify the `lb` and `ub` inputs
    !>  argument as `integer`s of default kind \IK, such that the random numbers are generated from the
    !>  processor-dependent character interval `[char(lb), char(ub)]`.<br>
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{getUnifRand}
    !>  \include{lineno} example/pm_distUnif/getUnifRand/main.F90
    !>  \compilef{getUnifRand}
    !>  \output{getUnifRand}
    !>  \include{lineno} example/pm_distUnif/getUnifRand/main.out.F90
    !>  \postproc{getUnifRand}
    !>  \include{lineno} example/pm_distUnif/getUnifRand/main.py
    !>  \vis{getUnifRand}
    !>  \image html pm_distUnif/getUnifRand/getUnifRand.CK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \todo
    !>  The current random `integer` generator uses a simple double precision `real` conversion to `integer` values.<br>
    !>  While this works fairly well for most use cases, it may biased for generating large random `integer` of kind \IK4.<br>
    !>  A future remedy should use Bitmask with Rejection as described [here](https://www.pcg-random.org/posts/bounded-rands.html).<br>
    !>  As of 2021, the use of double precision (64-bit) vs. single-precision for random number generation increases the
    !>  computational cost of the algorithms by about three times.<br>
    !>
    !>  \final{getUnifRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! RNGD

    interface getUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getUnifRandRNGDDD_D0_LK() result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDDD_D0_LK
#endif
        use pm_kind, only: LKG => LK
        logical(LKG)                                            :: rand
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_SK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
    end function
#endif

#if SK4_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_SK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
    end function
#endif

#if SK3_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_SK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
    end function
#endif

#if SK2_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_SK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
    end function
#endif

#if SK1_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_SK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_IK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
    end function
#endif

#if IK4_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_IK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
    end function
#endif

#if IK3_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_IK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
    end function
#endif

#if IK2_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_IK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
    end function
#endif

#if IK1_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_IK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_LK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
    end function
#endif

#if LK4_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_LK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
    end function
#endif

#if LK3_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_LK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
    end function
#endif

#if LK2_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_LK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
    end function
#endif

#if LK1_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_LK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_CK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
    end function
#endif

#if CK4_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_CK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
    end function
#endif

#if CK3_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_CK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
    end function
#endif

#if CK2_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_CK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
    end function
#endif

#if CK1_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_CK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_RK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_RK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_RK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_RK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getUnifRandRNGDLU_D0_RK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGDLU_D1_SK5(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGDLU_D1_SK4(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGDLU_D1_SK3(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGDLU_D1_SK2(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGDLU_D1_SK1(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGDLU_D1_IK5(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGDLU_D1_IK4(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGDLU_D1_IK3(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGDLU_D1_IK2(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGDLU_D1_IK1(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGDLU_D1_LK5(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGDLU_D1_LK4(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGDLU_D1_LK3(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGDLU_D1_LK2(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGDLU_D1_LK1(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGDLU_D1_CK5(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGDLU_D1_CK4(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGDLU_D1_CK3(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGDLU_D1_CK2(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGDLU_D1_CK1(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGDLU_D1_RK5(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGDLU_D1_RK4(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGDLU_D1_RK3(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGDLU_D1_RK2(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGDLU_D1_RK1(lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGDLU_D2_SK5(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGDLU_D2_SK4(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGDLU_D2_SK3(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGDLU_D2_SK2(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGDLU_D2_SK1(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGDLU_D2_IK5(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGDLU_D2_IK4(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGDLU_D2_IK3(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGDLU_D2_IK2(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGDLU_D2_IK1(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGDLU_D2_LK5(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGDLU_D2_LK4(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGDLU_D2_LK3(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGDLU_D2_LK2(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGDLU_D2_LK1(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGDLU_D2_CK5(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGDLU_D2_CK4(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGDLU_D2_CK3(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGDLU_D2_CK2(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGDLU_D2_CK1(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGDLU_D2_RK5(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGDLU_D2_RK4(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGDLU_D2_RK3(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGDLU_D2_RK2(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGDLU_D2_RK1(lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGDLU_D3_SK5(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGDLU_D3_SK4(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGDLU_D3_SK3(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGDLU_D3_SK2(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGDLU_D3_SK1(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGDLU_D3_IK5(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGDLU_D3_IK4(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGDLU_D3_IK3(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGDLU_D3_IK2(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGDLU_D3_IK1(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGDLU_D3_LK5(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGDLU_D3_LK4(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGDLU_D3_LK3(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGDLU_D3_LK2(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGDLU_D3_LK1(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGDLU_D3_CK5(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGDLU_D3_CK4(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGDLU_D3_CK3(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGDLU_D3_CK2(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGDLU_D3_CK1(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGDLU_D3_RK5(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGDLU_D3_RK4(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGDLU_D3_RK3(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGDLU_D3_RK2(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGDLU_D3_RK1(lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGDLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGF

    interface getUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getUnifRandRNGFDD_D0_LK(rng) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFDD_D0_LK
#endif
        use pm_kind, only: LKG => LK
        logical(LKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_SK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_SK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_SK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_SK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_SK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_IK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_IK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_IK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_IK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_IK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_LK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_LK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_LK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_LK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_LK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_CK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_CK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_CK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_CK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_CK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_RK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_RK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_RK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_RK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getUnifRandRNGFLU_D0_RK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGFLU_D1_SK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGFLU_D1_SK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGFLU_D1_SK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGFLU_D1_SK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGFLU_D1_SK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGFLU_D1_IK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGFLU_D1_IK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGFLU_D1_IK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGFLU_D1_IK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGFLU_D1_IK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGFLU_D1_LK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGFLU_D1_LK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGFLU_D1_LK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGFLU_D1_LK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGFLU_D1_LK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGFLU_D1_CK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGFLU_D1_CK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGFLU_D1_CK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGFLU_D1_CK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGFLU_D1_CK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGFLU_D1_RK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGFLU_D1_RK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGFLU_D1_RK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGFLU_D1_RK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGFLU_D1_RK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGFLU_D2_SK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGFLU_D2_SK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGFLU_D2_SK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGFLU_D2_SK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGFLU_D2_SK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGFLU_D2_IK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGFLU_D2_IK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGFLU_D2_IK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGFLU_D2_IK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGFLU_D2_IK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGFLU_D2_LK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGFLU_D2_LK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGFLU_D2_LK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGFLU_D2_LK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGFLU_D2_LK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGFLU_D2_CK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGFLU_D2_CK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGFLU_D2_CK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGFLU_D2_CK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGFLU_D2_CK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGFLU_D2_RK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGFLU_D2_RK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGFLU_D2_RK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGFLU_D2_RK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGFLU_D2_RK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGFLU_D3_SK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGFLU_D3_SK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGFLU_D3_SK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGFLU_D3_SK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGFLU_D3_SK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGFLU_D3_IK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGFLU_D3_IK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGFLU_D3_IK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGFLU_D3_IK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGFLU_D3_IK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGFLU_D3_LK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGFLU_D3_LK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGFLU_D3_LK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGFLU_D3_LK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGFLU_D3_LK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGFLU_D3_CK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGFLU_D3_CK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGFLU_D3_CK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGFLU_D3_CK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGFLU_D3_CK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGFLU_D3_RK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGFLU_D3_RK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGFLU_D3_RK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGFLU_D3_RK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGFLU_D3_RK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGFLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(rngf_type)         , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGS

    interface getUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getUnifRandRNGSDD_D0_LK(rng) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSDD_D0_LK
#endif
        use pm_kind, only: LKG => LK
        logical(LKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_SK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_SK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_SK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_SK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_SK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_IK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_IK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_IK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_IK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_IK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_LK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_LK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_LK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_LK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_LK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_CK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_CK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_CK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_CK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_CK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_RK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_RK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_RK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_RK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getUnifRandRNGSLU_D0_RK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGSLU_D1_SK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGSLU_D1_SK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGSLU_D1_SK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGSLU_D1_SK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGSLU_D1_SK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGSLU_D1_IK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGSLU_D1_IK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGSLU_D1_IK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGSLU_D1_IK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGSLU_D1_IK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGSLU_D1_LK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGSLU_D1_LK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGSLU_D1_LK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGSLU_D1_LK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGSLU_D1_LK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGSLU_D1_CK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGSLU_D1_CK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGSLU_D1_CK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGSLU_D1_CK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGSLU_D1_CK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGSLU_D1_RK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGSLU_D1_RK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGSLU_D1_RK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGSLU_D1_RK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGSLU_D1_RK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGSLU_D2_SK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGSLU_D2_SK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGSLU_D2_SK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGSLU_D2_SK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGSLU_D2_SK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGSLU_D2_IK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGSLU_D2_IK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGSLU_D2_IK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGSLU_D2_IK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGSLU_D2_IK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGSLU_D2_LK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGSLU_D2_LK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGSLU_D2_LK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGSLU_D2_LK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGSLU_D2_LK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGSLU_D2_CK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGSLU_D2_CK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGSLU_D2_CK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGSLU_D2_CK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGSLU_D2_CK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGSLU_D2_RK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGSLU_D2_RK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGSLU_D2_RK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGSLU_D2_RK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGSLU_D2_RK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGSLU_D3_SK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGSLU_D3_SK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGSLU_D3_SK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGSLU_D3_SK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGSLU_D3_SK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGSLU_D3_IK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGSLU_D3_IK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGSLU_D3_IK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGSLU_D3_IK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGSLU_D3_IK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGSLU_D3_LK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGSLU_D3_LK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGSLU_D3_LK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGSLU_D3_LK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGSLU_D3_LK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGSLU_D3_CK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGSLU_D3_CK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGSLU_D3_CK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGSLU_D3_CK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGSLU_D3_CK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGSLU_D3_RK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGSLU_D3_RK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGSLU_D3_RK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGSLU_D3_RK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGSLU_D3_RK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGSLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGX

    interface getUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getUnifRandRNGXDD_D0_LK(rng) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXDD_D0_LK
#endif
        use pm_kind, only: LKG => LK
        logical(LKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_SK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_SK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_SK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_SK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_SK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_IK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_IK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_IK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_IK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_IK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_LK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_LK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_LK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_LK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_LK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_CK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_CK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_CK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_CK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_CK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_RK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_RK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_RK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_RK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getUnifRandRNGXLU_D0_RK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGXLU_D1_SK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGXLU_D1_SK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGXLU_D1_SK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGXLU_D1_SK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGXLU_D1_SK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGXLU_D1_IK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGXLU_D1_IK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGXLU_D1_IK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGXLU_D1_IK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGXLU_D1_IK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGXLU_D1_LK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGXLU_D1_LK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGXLU_D1_LK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGXLU_D1_LK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGXLU_D1_LK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGXLU_D1_CK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGXLU_D1_CK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGXLU_D1_CK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGXLU_D1_CK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGXLU_D1_CK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGXLU_D1_RK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGXLU_D1_RK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGXLU_D1_RK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGXLU_D1_RK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGXLU_D1_RK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGXLU_D2_SK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGXLU_D2_SK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGXLU_D2_SK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGXLU_D2_SK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGXLU_D2_SK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGXLU_D2_IK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGXLU_D2_IK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGXLU_D2_IK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGXLU_D2_IK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGXLU_D2_IK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGXLU_D2_LK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGXLU_D2_LK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGXLU_D2_LK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGXLU_D2_LK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGXLU_D2_LK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGXLU_D2_CK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGXLU_D2_CK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGXLU_D2_CK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGXLU_D2_CK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGXLU_D2_CK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGXLU_D2_RK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGXLU_D2_RK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGXLU_D2_RK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGXLU_D2_RK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGXLU_D2_RK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGXLU_D3_SK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGXLU_D3_SK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGXLU_D3_SK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGXLU_D3_SK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGXLU_D3_SK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGXLU_D3_IK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGXLU_D3_IK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGXLU_D3_IK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGXLU_D3_IK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGXLU_D3_IK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGXLU_D3_LK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGXLU_D3_LK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGXLU_D3_LK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGXLU_D3_LK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGXLU_D3_LK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGXLU_D3_CK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGXLU_D3_CK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGXLU_D3_CK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGXLU_D3_CK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGXLU_D3_CK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGXLU_D3_RK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGXLU_D3_RK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGXLU_D3_RK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGXLU_D3_RK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGXLU_D3_RK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGXLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGG

    interface getUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getUnifRandRNGGDD_D0_LK(rng) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGDD_D0_LK
#endif
        use pm_kind, only: LKG => LK
        logical(LKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_SK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_SK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_SK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_SK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_SK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        character(max(len(lb,IK),len(ub,IK)),SKG)               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_IK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_IK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_IK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_IK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_IK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_LK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_LK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_LK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_LK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_LK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        logical(LKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_CK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_CK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_CK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_CK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_CK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        complex(CKG)                                            :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_RK5(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_RK4(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_RK3(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_RK2(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getUnifRandRNGGLU_D0_RK1(rng, lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        real(RKG)                                               :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGGLU_D1_SK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGGLU_D1_SK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGGLU_D1_SK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGGLU_D1_SK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGGLU_D1_SK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        character(len(lb,IK),SKG)                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGGLU_D1_IK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGGLU_D1_IK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGGLU_D1_IK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGGLU_D1_IK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGGLU_D1_IK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        integer(IKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGGLU_D1_LK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGGLU_D1_LK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGGLU_D1_LK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGGLU_D1_LK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGGLU_D1_LK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        logical(LKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGGLU_D1_CK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGGLU_D1_CK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGGLU_D1_CK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGGLU_D1_CK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGGLU_D1_CK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        complex(CKG)                                            :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGGLU_D1_RK5(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGGLU_D1_RK4(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGGLU_D1_RK3(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGGLU_D1_RK2(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGGLU_D1_RK1(rng, lb, ub, s1) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1
        real(RKG)                                               :: rand(s1)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGGLU_D2_SK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGGLU_D2_SK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGGLU_D2_SK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGGLU_D2_SK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGGLU_D2_SK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        character(len(lb,IK),SKG)                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGGLU_D2_IK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGGLU_D2_IK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGGLU_D2_IK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGGLU_D2_IK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGGLU_D2_IK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        integer(IKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGGLU_D2_LK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGGLU_D2_LK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGGLU_D2_LK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGGLU_D2_LK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGGLU_D2_LK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        logical(LKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGGLU_D2_CK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGGLU_D2_CK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGGLU_D2_CK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGGLU_D2_CK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGGLU_D2_CK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        complex(CKG)                                            :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGGLU_D2_RK5(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGGLU_D2_RK4(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGGLU_D2_RK3(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGGLU_D2_RK2(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGGLU_D2_RK1(rng, lb, ub, s1, s2) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2
        real(RKG)                                               :: rand(s1, s2)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getUnifRandRNGGLU_D3_SK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK4_ENABLED
    impure module function getUnifRandRNGGLU_D3_SK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK3_ENABLED
    impure module function getUnifRandRNGGLU_D3_SK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK2_ENABLED
    impure module function getUnifRandRNGGLU_D3_SK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if SK1_ENABLED
    impure module function getUnifRandRNGGLU_D3_SK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        character(len(lb,IK),SKG)                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure module function getUnifRandRNGGLU_D3_IK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK4_ENABLED
    impure module function getUnifRandRNGGLU_D3_IK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK3_ENABLED
    impure module function getUnifRandRNGGLU_D3_IK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK2_ENABLED
    impure module function getUnifRandRNGGLU_D3_IK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if IK1_ENABLED
    impure module function getUnifRandRNGGLU_D3_IK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        integer(IKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure module function getUnifRandRNGGLU_D3_LK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK4_ENABLED
    impure module function getUnifRandRNGGLU_D3_LK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK3_ENABLED
    impure module function getUnifRandRNGGLU_D3_LK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK2_ENABLED
    impure module function getUnifRandRNGGLU_D3_LK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if LK1_ENABLED
    impure module function getUnifRandRNGGLU_D3_LK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        logical(LKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getUnifRandRNGGLU_D3_CK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK4_ENABLED
    impure module function getUnifRandRNGGLU_D3_CK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK3_ENABLED
    impure module function getUnifRandRNGGLU_D3_CK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK2_ENABLED
    impure module function getUnifRandRNGGLU_D3_CK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if CK1_ENABLED
    impure module function getUnifRandRNGGLU_D3_CK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        complex(CKG)                                            :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRandRNGGLU_D3_RK5(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK4_ENABLED
    impure module function getUnifRandRNGGLU_D3_RK4(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK3_ENABLED
    impure module function getUnifRandRNGGLU_D3_RK3(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function

#endif

#if RK2_ENABLED
    impure module function getUnifRandRNGGLU_D3_RK2(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRandRNGGLU_D3_RK1(rng, lb, ub, s1, s2, s3) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRandRNGGLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                    :: lb, ub
        integer(IK)             , intent(in)                    :: s1, s2, s3
        real(RKG)                                               :: rand(s1, s2, s3)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a uniform random scalar or `contiguous` array of arbitrary rank of randomly uniformly distributed discrete
    !>  `logical`, `integer`, `character` value(s), or continuous `real` or `complex value(s) within the specified input range.
    !>
    !>  \param[inout]   rng     :   The input/output scalar of type,
    !>                              <ol>
    !>                                  <li>    [rngf_type](@ref pm_distUnif::rngf_type), or
    !>                                  <li>    [splitmix64_type](@ref pm_distUnif::splitmix64_type), or
    !>                                  <li>    [xoshiro256ssg_type](@ref pm_distUnif::xoshiro256ssg_type),
    !>                                  <li>    [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                              </ol>
    !>                              containing the user-specified random number generator algorithm to be used.<br>
    !>                              The user must initialize the object with the corresponding type constructors if non-deterministic RNG are desired.
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type))
    !>  \param[out]     rand    :   The output object of either<br>
    !>                              <ol>
    !>                                  <li>    scalar of type `character` of kind \SKALL of arbitrary `len` type parameter or, <br>
    !>                              </ol>
    !>                              or `contiguous` array of the rank, shape, and size as other array-like arguments, of either <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary `len` type parameter or, <br>
    !>                                  <li>    type `integer` of kind \IKALL or, <br>
    !>                                  <li>    type `complex` of kind \CKALL or, <br>
    !>                                  <li>    type `real` of kind \RKALL or, <br>
    !>                                  <li>    type `logical` of kind \LKALL, <br>
    !>                              </ol>
    !>                              containing the uniformly-distributed random output value.<br>
    !>                              <ol>
    !>                                  <li>    If `rand` is of type `logical`, its value is either `.false.` or `.true.`.<br>
    !>                                  <li>    If `rand` is of type `integer`, its value is in the interval `[lb, ub]`.<br>
    !>                                  <li>    If `rand` is of type `complex` or `real`, its value is in the interval `[lb, ub)`.<br>
    !>                                  <li>    If `rand` is of type `character`, its value is by default in the interval `[char(1), char(127)]`.<br>
    !>                              </ol>
    !>  \param[in]      lb      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `rand`, representing the lower bound of the Uniform distribution.<br>
    !>                              <ol>
    !>                                  <li>    If `rand` is of type `character`, then `len(rand) == len(lb)` must hold.<br>
    !>                                  <li>    If `rand` is of type `logical`, then `lb` must be `.false.` (there is no other possibility).<br>
    !>                              </ol>
    !>                              (**optional**, default = `char(1)`, `repeat(char(1), len(rand))`, `-huge(rand)`, `.false.`, `(0.,0.)`, or `0.`,
    !>                              for a scalar `character`, vector `character`, `integer`, `complex`, or `real` output `rand`, respectively.<br>
    !>                              It must be present <b>if and only if</b> `ub` is also present.)
    !>  \param[in]      ub      :   The input scalar (or array of the same shape as `rand`) of the same type and kind as `rand`,
    !>                              representing the upper bound of the Uniform distribution.<br>
    !>                              If `rand` is of type `character`, then `len(rand) == len(ub)` must hold.<br>
    !>                              If `rand` is of type `logical`, then `ub` does not exist as an input argument (it must not be present).<br>
    !>                              (**optional**, default = `char(127)`, `repeat(char(127), len(rand))`, `+huge(rand)`, `.true.`, `(1.,1.)`, or `1.`,
    !>                              for a scalar `character`, vector `character`, `integer`, `complex`, or `real` output `rand`, respectively.<br>
    !>                              It must be present <b>if and only if</b> `ub` is also present.)
    !>
    !>  \interface{setUnifRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnif, only: setUnifRand
    !>
    !>      call setUnifRand(rand)  ! `rand` can be any intrinsic type: character, complex, logical, integer, real.
    !>      call setUnifRand(rand, lb, ub) ! `rand` can be any intrinsic type: character, complex, logical, integer, real.
    !>
    !>      call setUnifRand(rng, rand) ! `rand` can be any intrinsic type: character, complex, logical, integer, real.
    !>      call setUnifRand(rng, rand, lb, ub) ! `rand` can be any intrinsic type: character, complex, logical, integer, real.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `len(lb) == len(ub) .or. len(lb) == 1 .or. len(ub) == 1` for the corresponding input arguments of type `character`.<br>
    !>  The condition `lb <= ub` must hold for the corresponding input arguments where
    !>  logical values are compared by the procedures of module [pm_logicalCompare](@ref pm_logicalCompare) and
    !>  complex values are compared by the procedures of module [pm_complexCompareAll](@ref pm_complexCompareAll).<br>
    !>  \vericons
    !>
    !>  \impure
    !>  The procedures of this generic interface are `pure` when the argument `rng` is present.<br>
    !>
    !>  \elemental
    !>  The procedures of this generic interface are non-`elemental` when the argument `rng` is present.<br>
    !>
    !>  \remark
    !>  The procedures under this generic interface are carefully designed to avoid possible overflow due to
    !>  specifying huge negative and positive `lb` and `ub` limits of type `integer`, `complex`, `real`.<br>
    !>  This is possible at the cost of making the random number generation slightly more expensive
    !>  (by a few CPU cycles, equivalent to and extra multiplication).<br>
    !>
    !>  \remark
    !>  It is expected that the condition `lb <= ub` if the two input arguments are specified by the user.<br>
    !>  However, this condition is neither enforced nor checked at runtime within the procedures.<br>
    !>
    !>  \note
    !>  By default random characters are generated from the ASCII collating sequence to ensure portability across compilers and platforms.<br>
    !>  If random uniform characters from the processor's collating sequence are desired, specify the `lb`
    !>  and `ub` inputs argument as `integer`s of default kind \IK, such that the random numbers are generated from the
    !>  processor-dependent character interval `[char(lb), char(ub)]`.<br>
    !>
    !>  \see
    !>  [rngf](@ref pm_distUnif::rngf)<br>
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getUnifCDF](@ref pm_distUnif::getUnifCDF)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifRandState](@ref pm_distUnif::getUnifRandState)<br>
    !>  [setUnifRandState](@ref pm_distUnif::setUnifRandState)<br>
    !>  [rngu_type](@ref pm_distUnif::rngu_type)<br>
    !>  [rngf_type](@ref pm_distUnif::rngf_type)<br>
    !>  [splitmix64_type](@ref pm_distUnif::splitmix64_type)<br>
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type)<br>
    !>  [getUnifRandStateSize](@ref pm_distUnif::getUnifRandStateSize)<br>
    !>
    !>  \example{setUnifRand}
    !>  \include{lineno} example/pm_distUnif/setUnifRand/main.F90
    !>  \compilef{setUnifRand}
    !>  \output{setUnifRand}
    !>  \include{lineno} example/pm_distUnif/setUnifRand/main.out.F90
    !>  \postproc{setUnifRand}
    !>  \include{lineno} example/pm_distUnif/setUnifRand/main.py
    !>  \vis{setUnifRand}
    !>  \image html pm_distUnif/setUnifRand/setUnifRand.IK.png width=700
    !>  \image html pm_distUnif/setUnifRand/setUnifRand.CK.png width=700
    !>  \image html pm_distUnif/setUnifRand/setUnifRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnif](@ref test_pm_distUnif)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3}
    !>  \desc
    !>  \gfortran yields an internal compiler error with the expression `rand = nint(temp, kind = IKG)`
    !>  in `pm_distUnif@routines@IK.inc.F90` file when `IKG => integer_kinds(5)` on WSL OS.<br>
    !>  \remedy
    !>  For now, the expression is replaced with `rand = int(0.5d0 + temp, kind = IKG)`.<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifx{2025.0.0 20241008}
    !>  \desc
    !>  \ifx{2025.0.0 20241008} cannot compile the following two lines of code in the include file `pm_distUnif@routines.inc.F90`.<br>
    !>  \code{.F90}
    !>
    !>      call setUnifRand(RNG rand%re, lb%re, ub%re)
    !>      call setUnifRand(RNG rand%im, lb%im, ub%im)
    !>
    !>  \endcode
    !>  Note that `ifort` can readily compile the above lines of code.<br>
    !>  Uncomment the above two lines to regenerate the compile-time error.<br>
    !>  \remedy
    !>  For now, these two lines are commented out for Intel compilers and replaced with the following.<br>
    !>  \code{.F90}
    !>
    !>      call setUnifRand(RNG rand%re, real(lb, CKG), real(ub, CKG))
    !>      call setUnifRand(RNG rand%im, aimag(lb), aimag(ub))
    !>
    !>  \endcode
    !>
    !>  \final{setUnifRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! RNGD

    interface setUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_SK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
    end subroutine
#endif

#if SK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_SK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
    end subroutine
#endif

#if SK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_SK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
    end subroutine
#endif

#if SK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_SK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
    end subroutine
#endif

#if SK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_SK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_IK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if IK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_IK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if IK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_IK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if IK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_IK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if IK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_IK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_LK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if LK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_LK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if LK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_LK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if LK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_LK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if LK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_LK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_CK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if CK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_CK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if CK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_CK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if CK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_CK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
    end subroutine
#endif

#if CK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_CK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDDD_D0_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_SK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
    end subroutine
#endif

#if SK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_SK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
    end subroutine
#endif

#if SK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_SK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
    end subroutine
#endif

#if SK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_SK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
    end subroutine
#endif

#if SK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_SK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_IK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if IK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_IK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if IK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_IK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if IK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_IK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if IK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_IK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_LK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if LK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_LK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if LK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_LK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if LK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_LK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if LK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_LK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_CK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if CK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_CK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if CK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_CK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if CK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_CK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

#if CK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_CK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_RK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_RK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_RK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_RK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setUnifRandRNGDLU_D0_RK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGDLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGF

    interface setUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFDD_D0_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if IK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if LK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if CK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setUnifRandRNGFLU_D0_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGFLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(rngf_type)         , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGS

    interface setUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSDD_D0_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGSLU_D0_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D1_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D2_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGSDD_D3_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGSLU_D3_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGSLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(splitmix64_type)   , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGX

    interface setUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGDD_D0_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGGLU_D0_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D1_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D2_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGGDD_D3_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGGLU_D3_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGGLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssg_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RNGS

    interface setUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXDD_D0_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setUnifRandRNGXLU_D0_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D1_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D2_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_SK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_SK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_SK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_SK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_SK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_IK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_IK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_IK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_IK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_IK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_LK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_LK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_LK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_LK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_LK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_CK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_CK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_CK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_CK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_CK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGXDD_D3_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_SK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_SK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_SK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_SK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_SK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(out)                   :: rand(:,:,:)
        character(*,SKG)        , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_IK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_IK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_IK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_IK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_IK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(out)                   :: rand(:,:,:)
        integer(IKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_LK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_LK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_LK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_LK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_LK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(out)                   :: rand(:,:,:)
        logical(LKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_CK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_CK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_CK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_CK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_CK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(out)                   :: rand(:,:,:)
        complex(CKG)            , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_RK5(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_RK4(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_RK3(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_RK2(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRandRNGXLU_D3_RK1(rng, rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRandRNGXLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:,:,:)
        real(RKG)               , intent(in)                    :: lb, ub
        type(xoshiro256ssw_type) , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distUnif
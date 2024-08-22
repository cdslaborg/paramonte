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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>GenExpGamma distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>GenExpGamma distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  A variable \f$X\f$ is said to be GenExpGamma-distributed if its PDF with **location** \f$-\infty < \log(\sigma) < +\infty\f$,
!>  **scale (inverse rate)** \f$\omega > 0\f$, and **shape** \f$\kappa > 0\f$ parameters is described by the following equation,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x | \kappa, \omega, \log(\sigma)) =
!>      \frac{1}{\omega \Gamma(\kappa)} ~
!>      \exp\Bigg( \kappa\bigg(\frac{x - \log(\sigma)}{\omega}\bigg) - \exp\bigg(\frac{x - \log(\sigma)}{\omega}\bigg) \Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>
!>  where \f$\eta = \frac{1}{\omega \Gamma(\kappa)}\f$ is the **normalization factor** of the PDF.<br>
!>  When \f$\sigma = 1\f$, the GenExpGamma PDF simplifies to the form,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x) =
!>      \frac{1}{\omega \Gamma(\kappa)} ~
!>      \exp\Bigg(\kappa \frac{x}{\omega} - \exp\bigg(\frac{x}{\omega}\bigg) \Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>
!>  If \f$(\sigma, \omega) = (1, 1)\f$, that is, \f$x\f$ is standardized, the GenExpGamma PDF further simplifies to the form,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x) =
!>      \frac{1}{\Gamma(\kappa)} ~
!>      \exp\Bigg(\kappa x - \exp(x) \Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>
!>  Setting the shape parameter to \f$\kappa = 1\f$ further simplifies the PDF to the form,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x) = \exp\Bigg(x - \exp(x)\Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>
!>  The parameter \f$\log(\sigma)\f$ determines the horizontal location of the mode of this unimodal GenExpGamma PDF.<br>
!>  While the PDF is always unimodal, the overall height and steepness of the PDF curve are determined by the values of \f$\kappa\f$ and \f$\omega\f$.<br>
!>  In addition, the tails of the PDF are **fat**, in the sense that the PDF decreases **algebraically** rather than decreasing exponentially for large values of \f$x\f$.<br>
!>  The term **GenExpGamma** is shorthand for **Exponential-Gamma**, and the GenExpGamma distribution is sometimes referred to as the **generalized extreme value distribution**.<br>
!>  The term also sometimes collectively refers to a number of qualitatively similar probability distributions.<br>
!>  The GenExpGamma distribution is mathematically defined to be the distribution of \f$\log(X)\f$ where \f$X\f$ follows a GenGamma distribution.<br>
!>  GenExpGamma distribution is also significant because of the statistical extreme value theorem,
!>  which states that the maximum of a sample of i.i.d. random variables can only converge in distribution to one of three possible distributions,<br>
!>  <ol>
!>      <li>    Gumbel distribution,
!>      <li>    Frechet distribution,
!>      <li>    Weibull distribution,
!>  </ol>
!>  all of which are special cases of GenExpGamma distribution.<br>
!>  The GenExpGamma distribution is used in various branches of science to model a number of phenomena,
!>  including wind speeds in meteorology and ocean current speeds in ocean engineering.<br>
!>
!>  The **CDF** of the GenExpGamma distribution over a strictly-positive support \f$x \in (0, +\infty)\f$ with the three (shape, shape, scale)
!>  parameters \f$(\kappa > 0, \omega > 0, \sigma > 0)\f$ is defined by the **regularized** Lower Incomplete Gamma function as,
!>  \f{eqnarray}{
!>      \large
!>      \mathrm{CDF}(x | \kappa, \sigma)
!>      & = & P\big(\kappa, \frac{x - \log(\sigma)}{\omega}\big) \\
!>      & = & \frac{1}{\Gamma(\kappa)} \int_0^{\exp(\frac{x - \log(\sigma)}{\omega})} ~ t^{\kappa - 1}{\mathrm e}^{-t} ~ dt ~,
!>  \f}
!>  where \f$\Gamma(\kappa)\f$ represents the Gamma function.<br>
!>
!>  \warning
!>  The **GenExpGamma** distribution is also frequently **mistakenly** called **GenLogGamma** , which is an entirely different distribution.<br>
!>
!>  \note
!>  The relationship between the [GenExpGamma](@ref pm_distGenExpGamma) and [GenGamma](@ref pm_distGenGamma) distributions
!>  is similar to that of the [Normal](@ref pm_distNorm) and [LogNormal](@ref pm_distLogNorm) distributions.<br>
!>  In other words, a better more consistent naming for the GenExpGamma and GenGamma distributions could have been *GenGamma* and *GenLogGamma* distributions,
!>  respectively, similar to [Normal](@ref pm_distNorm) and [LogNormal](@ref pm_distLogNorm) distributions.<br>
!>
!>  \see
!>  [pm_distGamma](@ref pm_distGamma)<br>
!>  [pm_distGenGamma](@ref pm_distGenGamma)<br>
!>  Wolfram Research (2010), GenExpGammaDistribution, Wolfram Language function, https://reference.wolfram.com/language/ref/GenExpGammaDistribution.html (updated 2016)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF, The runtime performance of [getGenExpGammaLogPDF](@ref pm_distGenExpGamma::getGenExpGammaLogPDF) vs. [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF)}
!>  \include{lineno} benchmark/pm_distGenExpGamma/getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF/main.F90
!>  \compilefb{getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF}
!>  \postprocb{getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF}
!>  \include{lineno} benchmark/pm_distGenExpGamma/getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF/main.py
!>  \visb{getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF}
!>  \image html benchmark/pm_distGenExpGamma/getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF/benchmark.getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF.runtime.png width=1000
!>  \image html benchmark/pm_distGenExpGamma/getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF/benchmark.getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF.runtime.ratio.png width=1000
!>  \moralb{getGenExpGammaLogPDF_vs_setGenExpGammaLogPDF}
!>      -#  The procedures under the generic interface [getGenExpGammaLogPDF](@ref pm_distGenExpGamma::getGenExpGammaLogPDF) are functions while
!>          the procedures under the generic interface [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than
!>          the subroutine interface when the input `x` size is small.<br>
!>          Otherwise, the difference appears to be marginal and insignificant in most practical situations.<br>
!>
!>  \test
!>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distGenExpGamma

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distGenExpGamma"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type GenExpGamma
    !>  as defined in the description of [pm_distGenExpGamma](@ref pm_distGenExpGamma).
    !>
    !>  \details
    !>  See the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma) for the definition of the GenExpGamma distribution.
    !>
    !>  \interface{distGenExpGamma_type}
    !>  \code{.F90}
    !>
    !>      use pm_distGenExpGamma, only: distGenExpGamma_type
    !>      type(distGenExpGamma_type) :: distGenExpGamma
    !>
    !>      distGenExpGamma = distGenExpGamma_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distGenExpGamma_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distGenExpGamma_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the Probability Density Function (PDF) of the GenExpGamma distribution.<br>
    !>
    !>  \details
    !>  The **normalization factor** of the GenExpGamma PDF is defined as \f$\eta = 1 / (\omega \Gamma(\kappa))\f$
    !>  where \f$(\omega, \kappa)\f$ are the scale and shape parameters of the GenExpGamma distribution, respectively.<br>
    !>  For more information, see the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma).<br>
    !>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              containing the inverse of the scale parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDFNF`              :   The output scalar or array of the same shape as array-like input arguments,
    !>                              of the same type, kind, and highest rank as the input arguments containing
    !>                              the natural logarithm of the normalization factor of the distribution.
    !>
    !>  \interface{getGenExpGammaLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenExpGamma, only: getGenExpGammaLogPDFNF
    !>
    !>      logPDFNF = getGenExpGammaLogPDFNF(kappa)
    !>      logPDFNF = getGenExpGammaLogPDFNF(kappa, invOmega)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the distribution.<br>
    !>  The logic behind pre-computing the normalization factor is to speed up the calculations since the normalization
    !>  factor contains `log()` and `log_gamma()` computations which are computationally costly operations.<br>
    !>  See the benchmark below for more details.<br>
    !>
    !>  \see
    !>  [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF)<br>
    !>  [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF)<br>
    !>
    !>  \example{getGenExpGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaLogPDFNF/main.F90
    !>  \compilef{getGenExpGammaLogPDFNF}
    !>  \output{getGenExpGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaLogPDFNF/main.out.F90
    !>  \postproc{getGenExpGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaLogPDFNF/main.py
    !>  \vis{getGenExpGammaLogPDFNF}
    !>  \image html pm_distGenExpGamma/getGenExpGammaLogPDFNF/getGenExpGammaLogPDFNF.RK.png width=700
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getGenExpGammaLogPDFNF, The runtime performance of [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF) with and without `scale`}
    !>  \include{lineno} benchmark/pm_distGenExpGamma/getGenExpGammaLogPDFNF/main.F90
    !>  \compilefb{getGenExpGammaLogPDFNF}
    !>  \postprocb{getGenExpGammaLogPDFNF}
    !>  \include{lineno} benchmark/pm_distGenExpGamma/getGenExpGammaLogPDFNF/main.py
    !>  \visb{getGenExpGammaLogPDFNF}
    !>  \image html benchmark/pm_distGenExpGamma/getGenExpGammaLogPDFNF/benchmark.getGenExpGammaLogPDFNF.runtime.png width=1000
    !>  \image html benchmark/pm_distGenExpGamma/getGenExpGammaLogPDFNF/benchmark.getGenExpGammaLogPDFNF.runtime.ratio.png width=1000
    !>  \moralb{getGenExpGammaLogPDFNF}
    !>      -#  The procedures under the generic interface [getGenExpGammaLogPDFNF](@ref pm_distGenExpGamma::getGenExpGammaLogPDFNF)
    !>          aim to precompute the natural logarithm of the normalization factor of the GenExpGamma PDF.<br>
    !>          Based on the results of this benchmark, passing this precomputed normalization factor to
    !>          [getGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF)
    !>          or [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF)
    !>          appears to lead to significant performance improvement (20-30 times)
    !>          when the GenExpGamma PDF procedures are to called repeatedly with the same
    !>          set of input distribution parameters (`kappa`, `sigma`).<br>
    !>          This speedup is entirely due to avoiding the redundant
    !>          repeated computations of `log_gamma()` and `log()` terms in
    !>          [getGenExpGammaLogPDFNF](@ref pm_distGenExpGamma::getGenExpGammaLogPDFNF).<br>
    !>      -#  The initial closing in the performance gap for small array sizes (of \f$\sim4\f$ or less elements)
    !>          is likely artificial and due to extremely small cost of `simpleAddition` procedure which is likely
    !>          less than the time resolution of the processor (\f$\sim10ns\f$ as of V2 of ParaMonte library).<br>
    !>
    !>  \test
    !>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGenExpGammaLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGenExpGammaLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKD_RK5(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKD_RK4(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKD_RK3(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKD_RK2(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKD_RK1(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKO_RK5(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKO_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKO_RK4(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKO_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKO_RK3(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKO_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKO_RK2(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKO_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenExpGammaLogPDFNFKO_RK1(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDFNFKO_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the GenExpGamma distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma) for more information on the PDF of the GenExpGamma distribution.<br>
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the values at which the log(PDF) must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the inverse of the scale parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  logSigma    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the location parameter (\f$\log(\sigma)\f$) of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as other array-valued input arguments
    !>                              of the same type and kind as `x` containing the PDF of the specified GenExpGamma distribution.
    !>
    !>  \interface{getGenExpGammaLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenExpGamma, only: getGenExpGammaLogPDF
    !>
    !>      logPDF = getGenExpGammaLogPDF(x, kappa = kappa, invOmega = invOmega, logSigma = logSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `kappa > 0` and `invOmega > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The nature of the definition of the GenExpGamma distribution as detailed in [pm_distGenExpGamma](@ref pm_distGenExpGamma)
    !>  requires the input absolute value of `x` (or its standardized value \f$y = \frac{x - \log(\sigma)}{\omega}\f$) to be reasonably
    !>  small (e.g., `x < 88.`, `x < 709.`, and `x < 11356` for the 32-bits, 64-bits, and 128-bits real `x` values) to avoid computational arithmetic overflow.<br>
    !>
    !>  \see
    !>  [setGenExpGammaLogPDF](@ref pm_distGenExpGamma::setGenExpGammaLogPDF)<br>
    !>
    !>  \example{getGenExpGammaLogPDF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaLogPDF/main.F90
    !>  \compilef{getGenExpGammaLogPDF}
    !>  \output{getGenExpGammaLogPDF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaLogPDF/main.out.F90
    !>  \postproc{getGenExpGammaLogPDF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaLogPDF/main.py
    !>  \vis{getGenExpGammaLogPDF}
    !>  \image html pm_distGenExpGamma/getGenExpGammaLogPDF/getGenExpGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGenExpGammaLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGenExpGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenExpGammaLogPDF_RK5(x, kappa, invOmega, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenExpGammaLogPDF_RK4(x, kappa, invOmega, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenExpGammaLogPDF_RK3(x, kappa, invOmega, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenExpGammaLogPDF_RK2(x, kappa, invOmega, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenExpGammaLogPDF_RK1(x, kappa, invOmega, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the GenExpGamma distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma) for more information on the PDF of the GenExpGamma distribution.<br>
    !>
    !>  \param[out] logPDF      :   The output scalar (or array of the same shape as other array-valued input arguments)
    !>                              of type `real` of kind \RKALL containing the PDF of the specified GenExpGamma distribution.<br>
    !>  \param[in]  x           :   The input scalar (or array of the same shape as other array-valued arguments), of the same type and kind as `logPDF`,
    !>                              containing the points at which the `log(PDF)` must be computed.<br>
    !>  \param[in]  logPDFNF    :   The input scalar (or array of the same shape as other array-valued arguments),
    !>                              containing the natural logarithm of the normalization factor (\f$\eta\f$) of the distribution.<br>
    !>                              The value of this argument can be obtained by calling [getGenExpGammaLogPDFNF](@ref pm_distGenExpGamma::getGenExpGammaLogPDFNF).<br>
    !>                              (**optional**, default = `getGenExpGammaLogPDFNF(kappa, invOmega)`. It must be present <b>if and only if</b> `kappa` is also present.)
    !>  \param[in]  kappa       :   The input scalar (or array of the same shape as other array-valued arguments), containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present if `invOmega` is also present.)
    !>  \param[in]  invOmega    :   The input scalar (or array of the same shape as other array-valued arguments), containing the inverse of the scale parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present if `logSigma` is also present.)
    !>  \param[in]  logSigma    :   The input scalar (or array of the same shape as other array-valued arguments), containing the location parameter (\f$\log(\sigma)\f$) of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>
    !>  \interface{setGenExpGammaLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenExpGamma, only: setGenExpGammaLogPDF
    !>
    !>      call setGenExpGammaLogPDF(logPDF, x)
    !>      call setGenExpGammaLogPDF(logPDF, x, logSigma)
    !>      call setGenExpGammaLogPDF(logPDF, x, logPDFNF, kappa, invOmega)
    !>      call setGenExpGammaLogPDF(logPDF, x, logPDFNF, kappa, invOmega, logSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `logPDFNF = getExpGammaLogPDFNF(kappa)`, `kappa > 0`, and `invOmega > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The nature of the definition of the GenExpGamma distribution as detailed in [pm_distGenExpGamma](@ref pm_distGenExpGamma)
    !>  requires the input absolute value of `x` (or its standardized value \f$y = \frac{x - \log(\sigma)}{\omega}\f$) to be reasonably
    !>  small (e.g., `x < 88.`, `x < 709.`, and `x < 11356` for the 32-bits, 64-bits, and 128-bits real `x` values) to avoid computational arithmetic overflow.<br>
    !>  Since \f$\lim_{x\rightarrow+\infty} \pi(x) \rightarrow +0\f$, the output \f$\log\big(\pi(x)\big)\f$ is set to `-huge(x)` to indicate \f$-\infty\f$, when `x` overflows the maximum value of the `real` kind.<br>
    !>  In summary, the overflow should not cause any problems in practice at runtime and is gracefully handled by the procedures under this generic interface.<br>
    !>
    !>  \see
    !>  [getGenExpGammaLogPDF](@ref pm_distGenExpGamma::getGenExpGammaLogPDF)<br>
    !>
    !>  \example{setGenExpGammaLogPDF}
    !>  \include{lineno} example/pm_distGenExpGamma/setGenExpGammaLogPDF/main.F90
    !>  \compilef{setGenExpGammaLogPDF}
    !>  \output{setGenExpGammaLogPDF}
    !>  \include{lineno} example/pm_distGenExpGamma/setGenExpGammaLogPDF/main.out.F90
    !>  \postproc{setGenExpGammaLogPDF}
    !>  \include{lineno} example/pm_distGenExpGamma/setGenExpGammaLogPDF/main.py
    !>  \vis{setGenExpGammaLogPDF}
    !>  \image html pm_distGenExpGamma/setGenExpGammaLogPDF/setGenExpGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setGenExpGammaLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGenExpGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFDDDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFDDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFDDDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFDDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFDDDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFDDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFDDDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFDDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFDDDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFDDDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKDD_RK5(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKDD_RK4(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKDD_RK3(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKDD_RK2(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKDD_RK1(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOD_RK5(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOD_RK4(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOD_RK3(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOD_RK2(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOD_RK1(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOS_RK5(logPDF, x, logPDFNF, kappa, invOmega, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, logSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOS_RK4(logPDF, x, logPDFNF, kappa, invOmega, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, logSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOS_RK3(logPDF, x, logPDFNF, kappa, invOmega, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, logSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOS_RK2(logPDF, x, logPDFNF, kappa, invOmega, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, logSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaLogPDFNKOS_RK1(logPDF, x, logPDFNF, kappa, invOmega, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaLogPDFNKOS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, logSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  GenExpGamma distribution for an input `x` within the support of the distribution \f$x \in (-\infty,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma) for more information on the GenExpGamma CDF.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the inverse scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  logSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>
    !>  \interface{getGenExpGammaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenExpGamma, only: getGenExpGammaCDF
    !>
    !>      cdf = getGenExpGammaCDF(x, kappa = kappa, invOmega = invOmega, logSigma = logSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `kappa > 0` and `invOmega > 0` must hold for the input values `(kappa, invOmega)`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The choice of `logSigma` vs. `sigma` as the input argument to the procedures of this generic interface may seem rather odd.<br>
    !>  This choice was made deliberately to increase the similarity and consistency of the functional interface here with the
    !>  subroutine interface [setGenExpGammaCDF](@ref pm_distGenExpGamma::setGenExpGammaCDF).<br>
    !>  Note that the primary reason for preferring `logSigma` vs. `sigma` as the choice of input argument is the enhanced computational efficiency of the routines.<br>
    !>  Additionally, `logSigma` matches the alternative **rate** interpretation of the second parameter of the GenExpGamma distribution.<br>
    !>
    !>  \see
    !>  [setGenExpGammaCDF](@ref pm_distGenExpGamma::setGenExpGammaCDF)<br>
    !>
    !>  \example{getGenExpGammaCDF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaCDF/main.F90
    !>  \compilef{getGenExpGammaCDF}
    !>  \output{getGenExpGammaCDF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaCDF/main.out.F90
    !>  \postproc{getGenExpGammaCDF}
    !>  \include{lineno} example/pm_distGenExpGamma/getGenExpGammaCDF/main.py
    !>  \vis{getGenExpGammaCDF}
    !>  \image html pm_distGenExpGamma/getGenExpGammaCDF/getGenExpGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGenExpGammaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGenExpGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenExpGammaCDF_RK5(x, kappa, invOmega, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenExpGammaCDF_RK4(x, kappa, invOmega, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenExpGammaCDF_RK3(x, kappa, invOmega, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenExpGammaCDF_RK2(x, kappa, invOmega, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenExpGammaCDF_RK1(x, kappa, invOmega, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenExpGammaCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, logSigma
        real(RKG)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  GenExpGamma distribution for an input `x` within the support of the distribution \f$x \in (-\infty,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma) for more information on the GenExpGamma CDF.
    !>
    !>  \param[out] cdf             :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                                  and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>  \param[in]  x               :   The input scalar or array of the same shape as other array like arguments,
    !>                                  of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa           :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It must be present if `invOmega` are also present.)
    !>  \param[in]  invOmega        :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the inverse scale parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It must be present if `logSigma` is also present.)
    !>  \param[in]  logSigma        :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the location parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It can be present only if all previous arguments are also present.)
    !>  \param[out] info            :   The output scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm.<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the CDF**.<br>
    !>                                  This is likely to happen if the input value for `kappa` is too large.<br>
    !>
    !>  \interface{setGenExpGammaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenExpGamma, only: setGenExpGammaCDF
    !>
    !>      call setGenExpGammaCDF(cdf, x, info)
    !>      call setGenExpGammaCDF(cdf, x, kappa, info)
    !>      call setGenExpGammaCDF(cdf, x, kappa, invOmega, info)
    !>      call setGenExpGammaCDF(cdf, x, kappa, invOmega, logSigma, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `kappa > 0` and `invOmega > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGenExpGammaCDF](@ref pm_distGenExpGamma::getGenExpGammaCDF)<br>
    !>
    !>  \example{setGenExpGammaCDF}
    !>  \include{lineno} example/pm_distGenExpGamma/setGenExpGammaCDF/main.F90
    !>  \compilef{setGenExpGammaCDF}
    !>  \output{setGenExpGammaCDF}
    !>  \include{lineno} example/pm_distGenExpGamma/setGenExpGammaCDF/main.out.F90
    !>  \postproc{setGenExpGammaCDF}
    !>  \include{lineno} example/pm_distGenExpGamma/setGenExpGammaCDF/main.py
    !>  \vis{setGenExpGammaCDF}
    !>  \image html pm_distGenExpGamma/setGenExpGammaCDF/setGenExpGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenExpGamma](@ref test_pm_distGenExpGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setGenExpGammaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGenExpGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFDDD_RK5(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFDDD_RK4(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFDDD_RK3(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFDDD_RK2(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFDDD_RK1(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFDDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFXKDD_RK5(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFXKDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFXKDD_RK4(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFXKDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFXKDD_RK3(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFXKDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFXKDD_RK2(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFXKDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFXKDD_RK1(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFXKDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOD_RK5(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOD_RK4(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOD_RK3(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOD_RK2(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOD_RK1(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOS_RK5(cdf, x, kappa, invOmega, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOS_RK4(cdf, x, kappa, invOmega, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOS_RK3(cdf, x, kappa, invOmega, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOS_RK2(cdf, x, kappa, invOmega, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenExpGammaCDFKOS_RK1(cdf, x, kappa, invOmega, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenExpGammaCDFKOS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distGenExpGamma
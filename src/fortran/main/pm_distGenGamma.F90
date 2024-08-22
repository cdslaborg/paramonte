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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>GenGamma distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>GenGamma distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  A variable \f$X\f$ is said to be **Generalized Gamma** (GenGamma) distributed if its PDF with the **scale** \f$0 < \sigma < +\infty\f$,
!>  **shape** \f$0 < \omega < +\infty\f$, and **shape** \f$0 < \kappa < +\infty\f$ parameters is described by the following equation,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x | \kappa, \omega, \sigma) =
!>      \frac{1}{\sigma \omega \Gamma(\kappa)} ~
!>      \bigg( \frac{x}{\sigma} \bigg)^{\frac{\kappa}{\omega} - 1} \exp\Bigg( -\bigg(\frac{x}{\sigma}\bigg)^{\frac{1}{\omega}} \Bigg) ~,~ 0 < x < \infty
!>  \f}
!>
!>  where \f$\eta = \frac{1}{\sigma \omega \Gamma(\kappa)}\f$ is the **normalization factor** of the PDF.<br>
!>  When \f$\sigma = 1\f$, the GenGamma PDF simplifies to the form,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x) =
!>      \frac{1}{\omega \Gamma(\kappa)} ~
!>      x^{\frac{\kappa}{\omega} - 1} \exp\Bigg( -x^{\frac{1}{\omega}} \Bigg) ~,~ 0 < x < \infty
!>  \f}
!>
!>  If \f$(\sigma, \omega) = (1, 1)\f$, the GenGamma PDF further simplifies to the form,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x) =
!>      \frac{1}{\Gamma(\kappa)} ~
!>      x^{\kappa - 1} \exp(-x) ~,~ 0 < x < \infty
!>  \f}
!>
!>  Setting the shape parameter to \f$\kappa = 1\f$ further simplifies the PDF to the Exponential distribution PDF with the scale parameter \f$\sigma = 1\f$,
!>
!>  \f{equation}{
!>      \large
!>      \pi(x) = \exp(x) ~,~ 0 < x < \infty
!>  \f}
!>
!>  <ol>
!>      <li>    The parameter \f$\sigma\f$ determines the scale of the GenGamma PDF.<br>
!>      <li>    When \f$\omega = 1\f$, the GenGamma PDF reduces to the PDF of the [Gamma distribution](@ref pm_distGamma).<br>
!>      <li>    When \f$\kappa = 1, \omega = 1\f$, the GenGamma PDF reduces to the PDF of the [Exponential distribution](@ref pm_distExp).<br>
!>  </ol>
!>
!>  The **CDF** of the Generalized Gamma distribution over a strictly-positive support \f$x \in (0, +\infty)\f$ with the three (shape, shape, scale)
!>  parameters \f$(\kappa > 0, \omega > 0, \sigma > 0)\f$ is defined by the **regularized** [Lower Incomplete Gamma function](@ref pm_mathGamma) as,
!>  \f{eqnarray}{
!>      \large
!>      \mathrm{CDF}(x | \kappa, \omega, \sigma)
!>      & = & P\bigg(\kappa, \big(\frac{x}{\sigma}\big)^{\frac{1}{\omega}} \bigg) \\
!>      & = & \frac{1}{\Gamma(\kappa)} \int_0^{\big(\frac{x}{\sigma}\big)^{\frac{1}{\omega}}} ~ t^{\kappa - 1}{\mathrm e}^{-t} ~ dt ~,
!>  \f}
!>
!>  where \f$\Gamma(\kappa)\f$ represents the Gamma function.<br>
!>
!>  The distribution mean is given by,
!>
!>  \f{equation}{
!>      \large
!>      \overline{x} = \frac{\Gamma\left(\kappa + \omega\right)}{\Gamma(\kappa)} \sigma ~.
!>  \f}
!>
!>  The distribution mode is given by,
!>
!>  \f{equation}{
!>      \large
!>      \widehat{x} =
!>      \begin{cases}
!>          \sigma \left( \kappa - \omega \right)^\omega ~~~ , ~~~ \omega < \kappa ~, \nonumber \\
!>          0 ~~~ , ~~~ \kappa \leq \omega ~.
!>      \end{cases}
!>  \f}
!>
!>  The distribution variance is given by,
!>
!>  \f{equation}{
!>      \large
!>      \mathrm{VAR}(x) = \sigma^2 \left[ \frac{\Gamma(\kappa + 2\omega)}{\Gamma(\kappa)} - \left( \frac{\Gamma(\kappa + \omega)}{\Gamma(\kappa)} \right)^2 \right] ~.
!>  \f}
!>
!>  \note
!>  The relationship between the [GenExpGamma](@ref pm_distGenExpGamma) and [GenGamma](@ref pm_distGenGamma) distributions
!>  is similar to that of the [Normal](@ref pm_distNorm) and [LogNormal](@ref pm_distLogNorm) distributions.<br>
!>  In other words, a better more consistent naming for the GenExpGamma and GenGamma distributions could have been *GenGamma* and *GenLogGamma* distributions,
!>  respectively, similar to [Normal](@ref pm_distNorm) and [LogNormal](@ref pm_distLogNorm) distributions.<br>
!>
!>  \see
!>  [pm_distGamma](@ref pm_distGamma)<br>
!>  [pm_distGenExpGamma](@ref pm_distGenExpGamma)<br>
!>  Stacy, E. W. (1962). A generalization of the gamma distribution. The Annals of mathematical statistics, 33(3), 1187-1192.<br>
!>  Wolfram Research (2010), GenGammaDistribution, Wolfram Language function, https://reference.wolfram.com/language/ref/GenGammaDistribution.html (updated 2016)<br>
!>
!>  \test
!>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distGenGamma

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distGenGamma"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type GenGamma
    !>  as defined in the description of [pm_distGenGamma](@ref pm_distGenGamma).
    !>
    !>  \details
    !>  See the documentation of [pm_distGenGamma](@ref pm_distGenGamma) for the definition of the GenGamma distribution.
    !>
    !>  \interface{distGenGamma_type}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: distGenGamma_type
    !>      type(distGenGamma_type) :: distGenGamma
    !>
    !>      distGenGamma = distGenGamma_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distGenGamma_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distGenGamma_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the Probability Density Function (PDF) of the GenGamma distribution.<br>
    !>
    !>  \details
    !>  The **normalization factor** of the GenGamma PDF is defined as \f$\eta = 1 / (\sigma \omega \Gamma(\kappa))\f$
    !>  where \f$(\omega, \kappa)\f$ are the scale and shape parameters of the GenGamma distribution, respectively.<br>
    !>  For more information, see the documentation of [pm_distGenGamma](@ref pm_distGenGamma).<br>
    !>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              containing the inverse of the second shape parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present if `invSigma` is also present.)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the inverse of the scale parameter (\f$\sigma\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDFNF`              :   The output scalar or array of the same shape as array-like input arguments,
    !>                              of the same type, kind, and highest rank as the input arguments containing
    !>                              the natural logarithm of the normalization factor of the distribution.
    !>
    !>  \interface{getGenGammaLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: getGenGammaLogPDFNF
    !>
    !>      logPDFNF = getGenGammaLogPDFNF(kappa)
    !>      logPDFNF = getGenGammaLogPDFNF(kappa, invOmega)
    !>      logPDFNF = getGenGammaLogPDFNF(kappa, invOmega, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `kappa > 0`, `invOmega > 0`, and `invSigma > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
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
    !>  [getGenGammaLogPDF](@ref pm_distGenGamma::getGenGammaLogPDF)<br>
    !>  [setGenGammaLogPDF](@ref pm_distGenGamma::setGenGammaLogPDF)<br>
    !>
    !>  \example{getGenGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaLogPDFNF/main.F90
    !>  \compilef{getGenGammaLogPDFNF}
    !>  \output{getGenGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaLogPDFNF/main.out.F90
    !>  \postproc{getGenGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaLogPDFNF/main.py
    !>  \vis{getGenGammaLogPDFNF}
    !>  \image html pm_distGenGamma/getGenGammaLogPDFNF/getGenGammaLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGenGammaLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGenGammaLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKDD_RK5(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKDD_RK4(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKDD_RK3(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKDD_RK2(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKDD_RK1(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOD_RK5(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOD_RK4(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOD_RK3(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOD_RK2(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOD_RK1(kappa, invOmega) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: kappa, invOmega
        real(RKG)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOS_RK5(kappa, invOmega, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOS_RK4(kappa, invOmega, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOS_RK3(kappa, invOmega, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOS_RK2(kappa, invOmega, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenGammaLogPDFNFKOS_RK1(kappa, invOmega, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDFNFKOS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the GenGamma distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenGamma](@ref pm_distGenGamma) for more information on the PDF of the GenGamma distribution.<br>
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the values at which the log(PDF) must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the inverse of the second shape parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the inverse of the scale parameter (\f$\sigma\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as other array-valued input arguments
    !>                              of the same type and kind as `x` containing the PDF of the specified GenGamma distribution.
    !>
    !>  \interface{getGenGammaLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: getGenGammaLogPDF
    !>
    !>      logPDF = getGenGammaLogPDF(x, kappa = kappa, invOmega = invOmega, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `kappa > 0`, `invOmega > 0`, and `invSigma > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGenGammaLogPDF](@ref pm_distGenGamma::setGenGammaLogPDF)<br>
    !>
    !>  \example{getGenGammaLogPDF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaLogPDF/main.F90
    !>  \compilef{getGenGammaLogPDF}
    !>  \output{getGenGammaLogPDF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaLogPDF/main.out.F90
    !>  \postproc{getGenGammaLogPDF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaLogPDF/main.py
    !>  \vis{getGenGammaLogPDF}
    !>  \image html pm_distGenGamma/getGenGammaLogPDF/getGenGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGenGammaLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGenGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenGammaLogPDF_RK5(x, kappa, invOmega, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenGammaLogPDF_RK4(x, kappa, invOmega, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenGammaLogPDF_RK3(x, kappa, invOmega, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenGammaLogPDF_RK2(x, kappa, invOmega, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenGammaLogPDF_RK1(x, kappa, invOmega, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, invOmega, invSigma
        real(RKG)                                   :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the GenGamma distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenGamma](@ref pm_distGenGamma) for more information on the PDF of the GenGamma distribution.<br>
    !>
    !>  \param[out] logPDF      :   The output scalar (or array of the same shape as other array-valued input arguments)
    !>                              of type `real` of kind \RKALL containing the PDF of the specified GenGamma distribution.<br>
    !>  \param[in]  x           :   The input scalar (or array of the same shape as other array-valued arguments), of the same type and kind as `logPDF`,
    !>                              containing the points at which the `log(PDF)` must be computed.<br>
    !>  \param[in]  logPDFNF    :   The input scalar (or array of the same shape as other array-valued arguments),
    !>                              containing the natural logarithm of the normalization factor (\f$\eta\f$) of the distribution.<br>
    !>                              See the documentation of [pm_distGenGamma](@ref pm_distGenGamma::getGenGammaLogPDFNF) for the definition and relevant benchmarks.<br>
    !>                              (**optional**, default = `0`. It must be present <b>if and only if</b> `kappa` is also present.)
    !>                              The value of this argument can be obtained by calling [getGenGammaLogPDFNF](@ref pm_distGenGamma::getGenGammaLogPDFNF).<br>
    !>  \param[in]  kappa       :   The input scalar (or array of the same shape as other array-valued arguments), containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present if `invOmega` is also present.)
    !>  \param[in]  invOmega    :   The input scalar (or array of the same shape as other array-valued arguments), containing the inverse of the second shape parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present if `invSigma` is also present.)
    !>  \param[in]  invSigma    :   The input scalar (or array of the same shape as other array-valued arguments), containing the inverse of the scale parameter (\f$\sigma\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \interface{setGenGammaLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: setGenGammaLogPDF
    !>
    !>      call setGenGammaLogPDF(logPDF, x)
    !>      call setGenGammaLogPDF(logPDF, x, logPDFNF, kappa)
    !>      call setGenGammaLogPDF(logPDF, x, logPDFNF, kappa, invOmega)
    !>      call setGenGammaLogPDF(logPDF, x, logPDFNF, kappa, invOmega, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `kappa > 0`, `invOmega > 0`, and `invSigma > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGenGammaLogPDF](@ref pm_distGenGamma::getGenGammaLogPDF)<br>
    !>
    !>  \example{setGenGammaLogPDF}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaLogPDF/main.F90
    !>  \compilef{setGenGammaLogPDF}
    !>  \output{setGenGammaLogPDF}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaLogPDF/main.out.F90
    !>  \postproc{setGenGammaLogPDF}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaLogPDF/main.py
    !>  \vis{setGenGammaLogPDF}
    !>  \image html pm_distGenGamma/setGenGammaLogPDF/setGenGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setGenGammaLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGenGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFDDDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFDDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFDDDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFDDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFDDDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFDDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFDDDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFDDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFDDDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFDDDD_RK1
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
    PURE elemental module subroutine setGenGammaLogPDFNKDD_RK5(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKDD_RK4(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKDD_RK3(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKDD_RK2(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKDD_RK1(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKDD_RK1
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
    PURE elemental module subroutine setGenGammaLogPDFNKOD_RK5(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOD_RK4(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOD_RK3(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOD_RK2(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOD_RK1(logPDF, x, logPDFNF, kappa, invOmega)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOD_RK1
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
    PURE elemental module subroutine setGenGammaLogPDFNKOS_RK5(logPDF, x, logPDFNF, kappa, invOmega, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOS_RK4(logPDF, x, logPDFNF, kappa, invOmega, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOS_RK3(logPDF, x, logPDFNF, kappa, invOmega, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOS_RK2(logPDF, x, logPDFNF, kappa, invOmega, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaLogPDFNKOS_RK1(logPDF, x, logPDFNF, kappa, invOmega, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaLogPDFNKOS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invOmega, invSigma
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
    !>  Generalized Gamma distribution for an input `x` within the support of the distribution \f$x \in (0,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenGamma](@ref pm_distGenGamma) for more information on the GenGamma CDF.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the inverse of the second shape parameter (\f$\omega\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>
    !>  \interface{getGenGammaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: getGenGammaCDF
    !>
    !>      cdf = getGenGammaCDF(x, kappa = kappa, invOmega = invOmega, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition \f$x \in (0,+\infty)\f$, `kappa > 0`, `invOmega > 0`, and `invSigma > 0` must hold for the input values `(x, kappa, invOmega, invSigma)`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGenGammaCDF](@ref pm_distGenGamma::setGenGammaCDF)<br>
    !>
    !>  \example{getGenGammaCDF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaCDF/main.F90
    !>  \compilef{getGenGammaCDF}
    !>  \output{getGenGammaCDF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaCDF/main.out.F90
    !>  \postproc{getGenGammaCDF}
    !>  \include{lineno} example/pm_distGenGamma/getGenGammaCDF/main.py
    !>  \vis{getGenGammaCDF}
    !>  \image html pm_distGenGamma/getGenGammaCDF/getGenGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGenGammaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGenGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGenGammaCDF_RK5(x, kappa, invOmega, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGenGammaCDF_RK4(x, kappa, invOmega, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGenGammaCDF_RK3(x, kappa, invOmega, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGenGammaCDF_RK2(x, kappa, invOmega, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGenGammaCDF_RK1(x, kappa, invOmega, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenGammaCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invOmega, invSigma
        real(RKG)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  Generalized Gamma distribution for an input `x` within the support of the distribution \f$x \in (0,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGenGamma](@ref pm_distGenGamma) for more information on the GenGamma CDF.
    !>
    !>  \param[out] cdf             :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                                  and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>  \param[in]  x               :   The input scalar or array of the same shape as other array like arguments,
    !>                                  of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa           :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It must be present if `invOmega` are also present.)
    !>  \param[in]  invOmega        :   The input scalar or array of the same shape as other array-valued arguments,
    !>                                  containing the inverse of the second shape parameter (\f$\omega\f$) of the distribution.<br>
    !>                                  (**optional**, default = `1.`)
    !>  \param[in]  invSigma        :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It can be present only if all previous arguments are also present.)
    !>  \param[out] info            :   The output scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm.<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the CDF**.<br>
    !>                                  This is likely to happen if the input value for `kappa` is too large.<br>
    !>
    !>  \interface{setGenGammaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: setGenGammaCDF
    !>
    !>      call setGenGammaCDF(cdf, x, info)
    !>      call setGenGammaCDF(cdf, x, kappa, info)
    !>      call setGenGammaCDF(cdf, x, kappa, invOmega, info)
    !>      call setGenGammaCDF(cdf, x, kappa, invOmega, invSigma, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition \f$x \in (0,+\infty)\f$ must hold for the corresponding input arguments.<br>
    !>  The condition `invSigma > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `invOmega > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `kappa > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGenGammaCDF](@ref pm_distGenGamma::getGenGammaCDF)<br>
    !>
    !>  \example{setGenGammaCDF}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaCDF/main.F90
    !>  \compilef{setGenGammaCDF}
    !>  \output{setGenGammaCDF}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaCDF/main.out.F90
    !>  \postproc{setGenGammaCDF}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaCDF/main.py
    !>  \vis{setGenGammaCDF}
    !>  \image html pm_distGenGamma/setGenGammaCDF/setGenGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setGenGammaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGenGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenGammaCDFDDD_RK5(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaCDFDDD_RK4(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaCDFDDD_RK3(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaCDFDDD_RK2(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaCDFDDD_RK1(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFDDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenGammaCDFKDD_RK5(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaCDFKDD_RK4(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaCDFKDD_RK3(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaCDFKDD_RK2(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaCDFKDD_RK1(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOD_RK5(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOD_RK4(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOD_RK3(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOD_RK2(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOD_RK1(cdf, x, kappa, invOmega, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOS_RK5(cdf, x, kappa, invOmega, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOS_RK4(cdf, x, kappa, invOmega, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOS_RK3(cdf, x, kappa, invOmega, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOS_RK2(cdf, x, kappa, invOmega, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGenGammaCDFKOS_RK1(cdf, x, kappa, invOmega, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaCDFKOS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invOmega, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar or array of arbitrary rank of GenGamma-distributed random values with the specified shape and scale
    !>  parameters \f$(\kappa, \omega, \sigma)\f$ of the Generalized Gamma distribution corresponding to the procedure arguments `(kappa, omega, sigma)`.
    !>
    !>  See the documentation of [pm_distGenGamma](@ref pm_distGenGamma) for more
    !>  information on the Probability Density Function (PDF) of the Generalized Gamma distribution.
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG for Gamma RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG for Gamma RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type), implying the use of the intrinsic Fortran URNG.)
    !>  \param[out]     rand    :   The output scalar or
    !>                              <ol>
    !>                                  <li>    array of rank `1`, or<br>
    !>                                  <li>    array of arbitrary rank if the `rng` argument is missing or set to [rngf_type](@ref pm_distUnif::rngf_type), or<br>
    !>                              </ol>
    !>                              of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL.
    !>                              </ol>
    !>                              On output, it contains GenGamma-distributed random value(s).<br>
    !>  \param[in]      kappa   :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `rand`,
    !>                              representing the \f$\kappa\f$ shape parameter of the Generalized Gamma distribution.<br>
    !>  \param[in]      sigma   :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `rand`,
    !>                              representing the \f$\omega\f$ shape parameter of the [Generalized Gamma distribution](@ref pm_distGenGamma).<br>
    !>  \param[in]      sigma   :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `rand`,
    !>                              representing the \f$\sigma\f$ scale parameter of the [Generalized Gamma distribution](@ref pm_distGenGamma).<br>
    !>
    !>  \interface{setGenGammaRand}
    !>  \code{.F90}
    !>
    !>      use pm_distGenGamma, only: setGenGammaRand
    !>
    !>      call setGenGammaRand(rand, kappa, omega, sigma)
    !>      call setGenGammaRand(rand(..), kappa, omega, sigma)
    !>      call setGenGammaRand(rng, rand, kappa, omega, sigma)
    !>      call setGenGammaRand(rng, rand(:), kappa, omega, sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < omega` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < sigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \recursive
    !>
    !>  \note
    !>  For repeated Gamma RNG with fixed `kappa`, it is best to pass a vector of `rand` to be filled
    !>  with random numbers rather than calling the procedures with scalar `rand` argument repeatedly.<br>
    !>  In addition to avoiding procedure call overhead, vectorized RGN in this particular case also avoids
    !>  an unnecessary division and square-root operation.<br>
    !>
    !>  \see
    !>  [getGenGammaLogPDF](@ref pm_distGenGamma::getGenGammaLogPDF)<br>
    !>  [setGenGammaLogPDF](@ref pm_distGenGamma::setGenGammaLogPDF)<br>
    !>  [getGenGammaCDF](@ref pm_distGenGamma::getGenGammaCDF)<br>
    !>  [setGenGammaCDF](@ref pm_distGenGamma::setGenGammaCDF)<br>
    !>
    !>  \example{setGenGammaRand}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaRand/main.F90
    !>  \compilef{setGenGammaRand}
    !>  \output{setGenGammaRand}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaRand/main.out.F90
    !>  \postproc{setGenGammaRand}
    !>  \include{lineno} example/pm_distGenGamma/setGenGammaRand/main.py
    !>  \vis{setGenGammaRand}
    !>  \image html pm_distGenGamma/setGenGammaRand/setGenGammaRand.RKG.png width=700
    !>
    !>  \test
    !>  [test_pm_distGenGamma](@ref test_pm_distGenGamma)
    !>
    !>  \final{setGenGammaRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGenGammaRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGenGammaRandRNGD_D0_RK5(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGenGammaRandRNGD_D0_RK4(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGenGammaRandRNGD_D0_RK3(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGenGammaRandRNGD_D0_RK2(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGenGammaRandRNGD_D0_RK1(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGenGammaRandRNGF_D0_RK5(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGenGammaRandRNGF_D0_RK4(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGenGammaRandRNGF_D0_RK3(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGenGammaRandRNGF_D0_RK2(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGenGammaRandRNGF_D0_RK1(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D0_RK5(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D0_RK4(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D0_RK3(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D0_RK2(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D0_RK1(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setGenGammaRandRNGD_D1_RK5(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setGenGammaRandRNGD_D1_RK4(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setGenGammaRandRNGD_D1_RK3(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setGenGammaRandRNGD_D1_RK2(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setGenGammaRandRNGD_D1_RK1(rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setGenGammaRandRNGF_D1_RK5(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setGenGammaRandRNGF_D1_RK4(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setGenGammaRandRNGF_D1_RK3(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setGenGammaRandRNGF_D1_RK2(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setGenGammaRandRNGF_D1_RK1(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGF_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D1_RK5(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D1_RK4(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D1_RK3(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D1_RK2(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGenGammaRandRNGX_D1_RK1(rng, rand, kappa, omega, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGenGammaRandRNGX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, omega, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distGenGamma
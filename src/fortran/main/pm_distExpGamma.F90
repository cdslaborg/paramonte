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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>ExpGamma distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>ExpGamma distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  A variable \f$X\f$ is said to be ExpGamma-distributed if its PDF with **location** \f$0 < \log(\sigma) < +\infty\f$,
!>  **scale (inverse rate)** \f$> 0\f$, and **shape** \f$\kappa > 0\f$ parameters is described by the following equation,
!>  \f{equation}{
!>      \large
!>      \pi(x | \kappa, \log(\sigma)) =
!>      \frac{1}{\Gamma(\kappa)} ~
!>      \exp\Bigg( \kappa\bigg(x - \log(\sigma)\bigg) - \exp\bigg(x - \log(\sigma)\bigg) \Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>  where \f$\eta = \frac{1}{\Gamma(\kappa)}\f$ is the **normalization factor** of the PDF.<br>
!>
!>  When \f$\sigma = 1\f$, the ExpGamma PDF simplifies to the form,
!>  \f{equation}{
!>      \large
!>      \pi(x) =
!>      \frac{1}{\Gamma(\kappa)} ~
!>      \exp\Bigg(\kappa x - \exp(x) \Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>
!>  Setting the shape parameter to \f$\kappa = 1\f$ further simplifies the PDF to the form,
!>  \f{equation}{
!>      \large
!>      \pi(x) = \exp\Bigg(x - \exp(x)\Bigg) ~,~ -\infty < x < \infty
!>  \f}
!>
!>  The parameter \f$\log(\sigma)\f$ determines the horizontal location of the mode of this unimodal ExpGamma PDF.<br>
!>  The ExpGamma distribution is mathematically defined to be the distribution of \f$\log(X)\f$ where \f$X\f$ follows a Gamma distribution.<br>
!>
!>  The **CDF** of the ExpGamma distribution over a strictly-positive support \f$x \in (0, +\infty)\f$ with the two (shape, scale)
!>  parameters \f$(\kappa > 0, \sigma > 0)\f$ is defined by the **regularized** Lower Incomplete Gamma function as,
!>  \f{eqnarray}{
!>      \large
!>      \mathrm{CDF}(x | \kappa, \sigma)
!>      & = & P\big(\kappa, x - \log(\sigma)\big) \\
!>      & = & \frac{1}{\Gamma(\kappa)} \int_0^{\exp(x - \log(\sigma))} ~ t^{\kappa - 1}{\mathrm e}^{-t} ~ dt ~,
!>  \f}
!>  where \f$\Gamma(\kappa)\f$ represents the Gamma function.<br>
!>
!>  \warning
!>  The **ExpGamma** distribution is also frequently **mistakenly** called **LogGamma** , which is an entirely different distribution.<br>
!>
!>  \note
!>  The relationship between the [ExpGamma](@ref pm_distExpGamma) and [Gamma](@ref pm_distGamma) distributions
!>  is similar to that of the [Normal](@ref pm_distNorm) and [LogNormal](@ref pm_distLogNorm) distributions.<br>
!>  In other words, a more consistent naming for the [ExpGamma](@ref pm_distExpGamma) and Gamma distributions could have been
!>  Gamma and LogGamma respectively, similar to [Normal](@ref pm_distNorm) and [LogNormal](@ref pm_distLogNorm) distributions.<br>
!>
!>  \see
!>  [pm_distGamma](@ref pm_distGamma)<br>
!>  [pm_distExpGamma](@ref pm_distExpGamma)<br>
!>  [pm_distGenGamma](@ref pm_distGenGamma)<br>
!>  [pm_distGenExpGamma](@ref pm_distGenExpGamma)<br>
!>
!>  \test
!>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distExpGamma

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distExpGamma"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type ExpGamma
    !>  as defined in the description of [pm_distExpGamma](@ref pm_distExpGamma).
    !>
    !>  \details
    !>  See the documentation of [pm_distExpGamma](@ref pm_distExpGamma) for the definition of the ExpGamma distribution.
    !>
    !>  \interface{distExpGamma_type}
    !>  \code{.F90}
    !>
    !>      use pm_distExpGamma, only: distExpGamma_type
    !>      type(distExpGamma_type) :: distExpGamma
    !>
    !>      distExpGamma = distExpGamma_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distExpGamma_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distExpGamma_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the Probability Density Function (PDF) of the ExpGamma distribution.<br>
    !>
    !>  \details
    !>  The **normalization factor** of the ExpGamma PDF is defined as \f$\eta = 1 / \Gamma(\kappa)\f$
    !>  where \f$\kappa\f$ is the shape parameter of the ExpGamma distribution.<br>
    !>  For more information, see the documentation of [pm_distExpGamma](@ref pm_distExpGamma).<br>
    !>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>
    !>  \return
    !>  `logPDFNF`              :   The output scalar or array of the same shape as array-like input arguments,
    !>                              of the same type, kind, and highest rank as the input arguments containing
    !>                              the natural logarithm of the normalization factor of the distribution.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExpGamma, only: getExpGammaLogPDFNF
    !>
    !>      logPDFNF = getExpGammaLogPDFNF(kappa)
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
    !>  [setExpGammaLogPDF](@ref pm_distExpGamma::setExpGammaLogPDF)<br>
    !>  [setExpGammaLogPDF](@ref pm_distExpGamma::setExpGammaLogPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaLogPDFNF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaLogPDFNF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaLogPDFNF/main.py
    !>  \vis
    !>  \image html pm_distExpGamma/getExpGammaLogPDFNF/getExpGammaLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getExpGammaLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpGammaLogPDFNF_RK5(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDFNF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpGammaLogPDFNF_RK4(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDFNF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpGammaLogPDFNF_RK3(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDFNF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpGammaLogPDFNF_RK2(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDFNF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpGammaLogPDFNF_RK1(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDFNF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: kappa
        real(RKG)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the ExpGamma distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distExpGamma](@ref pm_distExpGamma) for more information on the PDF of the ExpGamma distribution.<br>
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the values at which the log(PDF) must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  logSigma    :   The input scalar or array of the same shape as other array-valued arguments,
    !>                              containing the location parameter (\f$\log(\sigma)\f$) of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as other array-valued input arguments
    !>                              of the same type and kind as `x` containing the PDF of the specified ExpGamma distribution.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExpGamma, only: getExpGammaLogPDF
    !>
    !>      logPDF = getExpGammaLogPDF(x, kappa = kappa, logSigma = logSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `kappa > 0` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setExpGammaLogPDF](@ref pm_distExpGamma::setExpGammaLogPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaLogPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaLogPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaLogPDF/main.py
    !>  \vis
    !>  \image html pm_distExpGamma/getExpGammaLogPDF/getExpGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getExpGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpGammaLogPDF_RK5(x, kappa, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpGammaLogPDF_RK4(x, kappa, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpGammaLogPDF_RK3(x, kappa, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpGammaLogPDF_RK2(x, kappa, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpGammaLogPDF_RK1(x, kappa, logSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: x
        real(RKG)   , intent(in)    , optional      :: kappa, logSigma
        real(RKG)                                   :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the ExpGamma distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distExpGamma](@ref pm_distExpGamma) for more information on the PDF of the ExpGamma distribution.<br>
    !>
    !>  \param[out] logPDF      :   The output scalar (or array of the same shape as other array-valued input arguments)
    !>                              of type `real` of kind \RKALL containing the PDF of the specified ExpGamma distribution.<br>
    !>  \param[in]  x           :   The input scalar (or array of the same shape as other array-valued arguments), of the same type and kind as `logPDF`,
    !>                              containing the points at which the `log(PDF)` must be computed.<br>
    !>  \param[in]  logPDFNF    :   The input scalar (or array of the same shape as other array-valued arguments),
    !>                              containing the natural logarithm of the normalization factor (\f$\eta\f$) of the distribution.<br>
    !>                              The value of this argument can be obtained by calling [getExpGammaLogPDFNF](@ref pm_distExpGamma::getExpGammaLogPDFNF).<br>
    !>                              (**optional**, default = `getExpGammaLogPDFNF(kappa)`. It must be present <b>if and only if</b> `kappa` is also present.)
    !>  \param[in]  kappa       :   The input scalar (or array of the same shape as other array-valued arguments), containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present if `logSigma` is also present.)
    !>  \param[in]  logSigma    :   The input scalar (or array of the same shape as other array-valued arguments), containing the location parameter (\f$\log(\sigma)\f$) of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExpGamma, only: setExpGammaLogPDF
    !>
    !>      call setExpGammaLogPDF(logPDF, x)
    !>      call setExpGammaLogPDF(logPDF, x, logPDFNF, kappa)
    !>      call setExpGammaLogPDF(logPDF, x, logPDFNF, kappa, logSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `logPDFNF = getExpGammaLogPDFNF(kappa)` and `kappa > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setExpGammaLogPDF](@ref pm_distExpGamma::setExpGammaLogPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExpGamma/setExpGammaLogPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExpGamma/setExpGammaLogPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExpGamma/setExpGammaLogPDF/main.py
    !>  \vis
    !>  \image html pm_distExpGamma/setExpGammaLogPDF/setExpGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setExpGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFDDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFDDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFDDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFDDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFDDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFDDD_RK1
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
    PURE elemental module subroutine setExpGammaLogPDFNKD_RK5(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKD_RK4(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKD_RK3(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKD_RK2(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKD_RK1(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKD_RK1
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
    PURE elemental module subroutine setExpGammaLogPDFNKS_RK5(logPDF, x, logPDFNF, kappa, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, logSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKS_RK4(logPDF, x, logPDFNF, kappa, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, logSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKS_RK3(logPDF, x, logPDFNF, kappa, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, logSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKS_RK2(logPDF, x, logPDFNF, kappa, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, logSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpGammaLogPDFNKS_RK1(logPDF, x, logPDFNF, kappa, logSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaLogPDFNKS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, logSigma
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
    !>  ExpGamma distribution for an input `x` within the support of the distribution \f$x \in (-\infty,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distExpGamma](@ref pm_distExpGamma) for more information on the ExpGamma CDF.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  logSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExpGamma, only: getExpGammaCDF
    !>
    !>      cdf = getExpGammaCDF(x, kappa = kappa, logSigma = logSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `kappa > 0` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The choice of `logSigma` vs. `sigma` as the input argument to the procedures of this generic interface may seem rather odd.<br>
    !>  This choice was made deliberately to increase the similarity and consistency of the functional interface here with the
    !>  subroutine interface [setExpGammaCDF](@ref pm_distExpGamma::setExpGammaCDF).<br>
    !>  Note that the primary reason for preferring `logSigma` vs. `sigma` as the choice of input argument is the enhanced computational efficiency of the routines.<br>
    !>  Additionally, `logSigma` matches the alternative **rate** interpretation of the second parameter of the ExpGamma distribution.<br>
    !>
    !>  \see
    !>  [setExpGammaCDF](@ref pm_distExpGamma::setExpGammaCDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaCDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaCDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExpGamma/getExpGammaCDF/main.py
    !>  \vis
    !>  \image html pm_distExpGamma/getExpGammaCDF/getExpGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getExpGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpGammaCDF_RK5(x, kappa, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpGammaCDF_RK4(x, kappa, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpGammaCDF_RK3(x, kappa, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpGammaCDF_RK2(x, kappa, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, logSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpGammaCDF_RK1(x, kappa, logSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpGammaCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, logSigma
        real(RKG)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  ExpGamma distribution for an input `x` within the support of the distribution \f$x \in (-\infty,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distExpGamma](@ref pm_distExpGamma) for more information on the ExpGamma CDF.
    !>
    !>  \param[out] cdf             :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                                  and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>  \param[in]  x               :   The input scalar or array of the same shape as other array like arguments,
    !>                                  of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa           :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
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
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExpGamma, only: setExpGammaCDF
    !>
    !>      call setExpGammaCDF(cdf, x, info)
    !>      call setExpGammaCDF(cdf, x, kappa, info)
    !>      call setExpGammaCDF(cdf, x, kappa, logSigma, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `kappa > 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getExpGammaCDF](@ref pm_distExpGamma::getExpGammaCDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExpGamma/setExpGammaCDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExpGamma/setExpGammaCDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExpGamma/setExpGammaCDF/main.py
    !>  \vis
    !>  \image html pm_distExpGamma/setExpGammaCDF/setExpGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExpGamma](@ref test_pm_distExpGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setExpGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpGammaCDFDD_RK5(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpGammaCDFDD_RK4(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpGammaCDFDD_RK3(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpGammaCDFDD_RK2(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpGammaCDFDD_RK1(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpGammaCDFKD_RK5(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpGammaCDFKD_RK4(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpGammaCDFKD_RK3(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpGammaCDFKD_RK2(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpGammaCDFKD_RK1(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpGammaCDFKS_RK5(cdf, x, kappa, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpGammaCDFKS_RK4(cdf, x, kappa, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpGammaCDFKS_RK3(cdf, x, kappa, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpGammaCDFKS_RK2(cdf, x, kappa, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpGammaCDFKS_RK1(cdf, x, kappa, logSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpGammaCDFKS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, logSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distExpGamma
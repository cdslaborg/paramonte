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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Gamma distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Gamma distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the Gamma distribution over a strictly-positive support \f$x \in (0, +\infty)\f$
!>  is defined with the two (shape, scale) parameters \f$(\kappa > 0, \sigma > 0)\f$ as,
!>  \f{equation}{
!>
!>      \large
!>      \pi(x | \kappa, \sigma) = \frac {1} {\sigma\Gamma(\kappa)} \bigg(\frac{x}{\sigma}\bigg)^{\kappa - 1} \exp\bigg(-\frac{x}{\sigma}\bigg) ~,
!>
!>  \f}
!>
!>  where \f$\eta = \frac{1}{\sigma\Gamma(\kappa)}\f$ is the **normalization factor** of the distribution with \f$\Gamma(\kappa)\f$
!>  representing the Gamma function whose natural logarithm is returned by the Fortran intrinsic `log_gamma()`.<br>
!>  Note that \f$\lim_{x\to0} \pi(x | 0 < \kappa < 1, \sigma) \to +\infty\f$ \f$\f$.<br>
!>  However, this divergence is integrable and the Gamma PDF is properly normalized.<br>
!>
!>  The **CDF** of the Gamma distribution over a strictly-positive support \f$x \in (0, +\infty)\f$ with the two (shape, scale)
!>  parameters \f$(\kappa > 0, \sigma > 0)\f$ is defined by the **regularized** Lower Incomplete Gamma function as,
!>  \f{eqnarray}{
!>      \large
!>      \mathrm{CDF}(x | \kappa, \sigma)
!>      & = & P(\kappa, \frac{x}{\sigma}) \\
!>      & = & \frac{1}{\Gamma(\kappa)} \int_0^{\frac{x}{\sigma}} ~ t^{\kappa - 1}{\mathrm e}^{-t} ~ dt ~,
!>  \f}
!>  where \f$\Gamma(\kappa)\f$ represents the Gamma function.<br>
!>
!>  \note
!>  The **Cumulative Distribution Function (CDF)** of the Gamma distribution for a triple \f$(x, \kappa, \sigma)\f$ is
!>  directly returned by calling [getGammaIncLow(sigma * x, kappa)](@ref pm_mathGamma::getGammaIncLow).<br>
!>
!>  \see
!>  Press et al., 1992, Numerical Recipes.<br>
!>  [Marsaglia and Tsang, 2000, A simple method for generating gamma variables](https://dl.acm.org/doi/10.1145/358407.358414)<br>
!>  [Intel Math Kernel Library](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top.html)<br>
!>  [Blog Post by Sukhbinder Singh](https://sukhbinder.wordpress.com/fortran-random-number-generation/)<br>
!>  [Blog Post by Yutaka Masuda](https://masuday.github.io/fortran_tutorial/random.html)<br>
!>  [cdflib](https://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html)<br>
!>
!>  \test
!>  [test_pm_distGamma](@ref test_pm_distGamma)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distGamma

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distGamma"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Gamma
    !>  as defined in the description of [pm_distGamma](@ref pm_distGamma).
    !>
    !>  \details
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for the definition of the Gamma distribution.
    !>
    !>  \interface{distGamma_type}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: distGamma_type
    !>      type(distGamma_type) :: distGamma
    !>
    !>      distGamma = distGamma_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distGamma_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distGamma_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the Probability Density Function (PDF)
    !>  of the Gamma distribution for an input parameter set \f$(\kappa,\sigma)\f$.
    !>
    !>  \brief
    !>  The natural logarithm of the normalization factor of the Gamma distribution is given by,
    !>  \f{equation}{
    !>
    !>      \large
    !>      \ms{logPDFNF} \equiv \log(\eta) = \log\bigg( \frac {1} {\sigma\Gamma(\kappa)} \bigg) ~,
    !>
    !>  \f}
    !>
    !>  where \f$\Gamma(\kappa)\f$ is the Gamma function whose natural logarithm is returned by the Fortran intrinsic `log_gamma()`.<br>
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for more information on the Gamma distribution.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factor of the Gamma distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the Gamma PDF to improve the runtime performance by eliminating redundant calculations.
    !>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of type `real` of kind \RKALL, containing the shape parameter of the distribution.<br>
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `kappa`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDFNF`              :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `kappa`, containing the natural logarithm
    !>                              of the normalization factor of the PDF of Gamma distribution.<br>
    !>
    !>  \interface{getGammaLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: getGammaLogPDFNF
    !>
    !>      logPDFNF = getGammaLogPDFNF(kappa)
    !>      logPDFNF = getGammaLogPDFNF(kappa, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGammaLogPDF](@ref pm_distGamma::getGammaLogPDF)<br>
    !>  [setGammaLogPDF](@ref pm_distGamma::setGammaLogPDF)<br>
    !>
    !>  \example{getGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGamma/getGammaLogPDFNF/main.F90
    !>  \compilef{getGammaLogPDFNF}
    !>  \output{getGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGamma/getGammaLogPDFNF/main.out.F90
    !>  \postproc{getGammaLogPDFNF}
    !>  \include{lineno} example/pm_distGamma/getGammaLogPDFNF/main.py
    !>  \vis{getGammaLogPDFNF}
    !>  \image html pm_distGamma/getGammaLogPDFNF/getGammaLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGammaLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGammaLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGammaLogPDFNFKD_RK5(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: kappa
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGammaLogPDFNFKD_RK4(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: kappa
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGammaLogPDFNFKD_RK3(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: kappa
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGammaLogPDFNFKD_RK2(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: kappa
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGammaLogPDFNFKD_RK1(kappa) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: kappa
        real(RKG)                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGammaLogPDFNFKS_RK5(kappa, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: kappa, invSigma
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGammaLogPDFNFKS_RK4(kappa, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: kappa, invSigma
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGammaLogPDFNFKS_RK3(kappa, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: kappa, invSigma
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGammaLogPDFNFKS_RK2(kappa, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: kappa, invSigma
        real(RKG)                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGammaLogPDFNFKS_RK1(kappa, invSigma) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDFNFKS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: kappa, invSigma
        real(RKG)                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the
    !>  Gamma distribution for an input `x` within the support of the distribution \f$x \in (0,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for more information on the Gamma distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the PDF must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.<br>
    !>
    !>  \interface{getGammaLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: getGammaLogPDF
    !>
    !>      logPDF = getGammaLogPDF(x, kappa = kappa, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The choice of `invSigma` vs. `sigma` as the input argument to the procedures of this generic interface may seem rather odd.<br>
    !>  This choice was made deliberately to increase the similarity and consistency of the functional interface here with the
    !>  subroutine interface [setGammaLogPDF](@ref pm_distGamma::setGammaLogPDF).<br>
    !>  Note that the primary reason for preferring `invSigma` vs. `sigma` as the choice of input argument is the enhanced computational efficiency of the routines.<br>
    !>  Additionally, `invSigma` matches the alternative **rate** interpretation of the second parameter of the Gamma PDF.<br>
    !>
    !>  \note
    !>  See [setGammaLogPDF](@ref pm_distGamma::setGammaLogPDF) for a less flexible but more performant interface to the same functionality of this interface.
    !>
    !>  \see
    !>  [setGammaLogPDF](@ref pm_distGamma::setGammaLogPDF)<br>
    !>
    !>  \example{getGammaLogPDF}
    !>  \include{lineno} example/pm_distGamma/getGammaLogPDF/main.F90
    !>  \compilef{getGammaLogPDF}
    !>  \output{getGammaLogPDF}
    !>  \include{lineno} example/pm_distGamma/getGammaLogPDF/main.out.F90
    !>  \postproc{getGammaLogPDF}
    !>  \include{lineno} example/pm_distGamma/getGammaLogPDF/main.py
    !>  \vis{getGammaLogPDF}
    !>  \image html pm_distGamma/getGammaLogPDF/getGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGammaLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGammaLogPDF_RK5(x, kappa, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGammaLogPDF_RK4(x, kappa, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGammaLogPDF_RK3(x, kappa, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGammaLogPDF_RK2(x, kappa, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGammaLogPDF_RK1(x, kappa, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the
    !>  Gamma distribution for an input `x` within the support of the distribution \f$x \in (0,+\infty)\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for more information on the Gamma distribution.<br>
    !>
    !>  \param[out] logPDF      :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.<br>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the PDF must be computed.<br>
    !>  \param[in]  logPDFNF    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the natural logarithm of the normalization factor of the Gamma distribution.<br>
    !>                              The value of this input argument is readily returned by [getGammaLogPDFNF(kappa, invSigma)](@ref pm_distGamma::getGammaLogPDFNF).<br>
    !>                              significantly improves the runtime performance.<br>
    !>                              (**optional**, default = `0`. It must be present <b>if and only if</b> `kappa` is also present.)
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present <b>if and only if</b> `logPDFNF` is also present.)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`. It can be present only if all previous arguments are also present.)
    !>
    !>  \interface{setGammaLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: setGammaLogPDF
    !>
    !>      call setGammaLogPDF(logPDF, x)
    !>      call setGammaLogPDF(logPDF, x, logPDFNF, kappa)
    !>      call setGammaLogPDF(logPDF, x, logPDFNF, kappa, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  The condition `logPDFNF = getGammaLogPDFNF(kappa, invSigma)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  See [getGammaLogPDF](@ref pm_distGamma::getGammaLogPDF) for a more flexible but less performant interface to the same functionality of this interface.<br>
    !>
    !>  \see
    !>  [setGammaLogPDF](@ref pm_distGamma::setGammaLogPDF)<br>
    !>
    !>  \example{setGammaLogPDF}
    !>  \include{lineno} example/pm_distGamma/setGammaLogPDF/main.F90
    !>  \compilef{setGammaLogPDF}
    !>  \output{setGammaLogPDF}
    !>  \include{lineno} example/pm_distGamma/setGammaLogPDF/main.out.F90
    !>  \postproc{setGammaLogPDF}
    !>  \include{lineno} example/pm_distGamma/setGammaLogPDF/main.py
    !>  \vis{setGammaLogPDF}
    !>  \image html pm_distGamma/setGammaLogPDF/setGammaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setGammaLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGammaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaLogPDFDDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaLogPDFDDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaLogPDFDDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaLogPDFDDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaLogPDFDDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFDDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKD_RK5(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKD_RK4(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKD_RK3(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKD_RK2(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKD_RK1(logPDF, x, logPDFNF, kappa)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKS_RK5(logPDF, x, logPDFNF, kappa, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKS_RK4(logPDF, x, logPDFNF, kappa, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKS_RK3(logPDF, x, logPDFNF, kappa, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKS_RK2(logPDF, x, logPDFNF, kappa, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaLogPDFNKS_RK1(logPDF, x, logPDFNF, kappa, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaLogPDFNKS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, logPDFNF, kappa, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  Gamma distribution for an input `x` within the support of the distribution \f$x \in (0,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for more information on the Gamma CDF.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>
    !>  \interface{getGammaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: getGammaCDF
    !>
    !>      cdf = getGammaCDF(x, kappa = kappa, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The choice of `invSigma` vs. `sigma` as the input argument to the procedures of this generic interface may seem rather odd.<br>
    !>  This choice was made deliberately to increase the similarity and consistency of the functional interface here with the
    !>  subroutine interface [setGammaCDF](@ref pm_distGamma::setGammaCDF).<br>
    !>  Note that the primary reason for preferring `invSigma` vs. `sigma` as the choice of input argument is the enhanced computational efficiency of the routines.<br>
    !>  Additionally, `invSigma` matches the alternative **rate** interpretation of the second parameter of the Gamma distribution.<br>
    !>
    !>  \see
    !>  [setGammaCDF](@ref pm_distGamma::setGammaCDF)<br>
    !>
    !>  \example{getGammaCDF}
    !>  \include{lineno} example/pm_distGamma/getGammaCDF/main.F90
    !>  \compilef{getGammaCDF}
    !>  \output{getGammaCDF}
    !>  \include{lineno} example/pm_distGamma/getGammaCDF/main.out.F90
    !>  \postproc{getGammaCDF}
    !>  \include{lineno} example/pm_distGamma/getGammaCDF/main.py
    !>  \vis{getGammaCDF}
    !>  \image html pm_distGamma/getGammaCDF/getGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getGammaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGammaCDF_RK5(x, kappa, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGammaCDF_RK4(x, kappa, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGammaCDF_RK3(x, kappa, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGammaCDF_RK2(x, kappa, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGammaCDF_RK1(x, kappa, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: kappa, invSigma
        real(RKG)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  Gamma distribution for an input `x` within the support of the distribution \f$x \in (0,+\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for more information on the Gamma CDF.
    !>
    !>  \param[out] cdf             :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                                  and kind the input argument `x`, containing the CDF of the distribution at the specified `x`.<br>
    !>  \param[in]  x               :   The input scalar or array of the same shape as other array like arguments,
    !>                                  of type `real` of kind \RKALL, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  kappa           :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the shape parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It must be present if `invSigma` is also present.)
    !>  \param[in]  invSigma        :   The input scalar or array of the same shape as other array-like arguments,
    !>                                  of the same type and kind as `x`, containing the rate (inverse scale) parameter of the distribution.<br>
    !>                                  (**optional**, default = `1.`. It can be present only if all previous arguments are also present.)
    !>  \param[out] info            :   The output scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm.<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the CDF**.<br>
    !>                                  This is likely to happen if the input value for `kappa` is too large.<br>
    !>
    !>  \interface{setGammaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: setGammaCDF
    !>
    !>      call setGammaCDF(cdf, x, info)
    !>      call setGammaCDF(cdf, x, kappa, invSigma, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGammaCDF](@ref pm_distGamma::getGammaCDF)<br>
    !>
    !>  \example{setGammaCDF}
    !>  \include{lineno} example/pm_distGamma/setGammaCDF/main.F90
    !>  \compilef{setGammaCDF}
    !>  \output{setGammaCDF}
    !>  \include{lineno} example/pm_distGamma/setGammaCDF/main.out.F90
    !>  \postproc{setGammaCDF}
    !>  \include{lineno} example/pm_distGamma/setGammaCDF/main.py
    !>  \vis{setGammaCDF}
    !>  \image html pm_distGamma/setGammaCDF/setGammaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setGammaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGammaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaCDFDD_RK5(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaCDFDD_RK4(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaCDFDD_RK3(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaCDFDD_RK2(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaCDFDD_RK1(cdf, x, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaCDFKD_RK5(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaCDFKD_RK4(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaCDFKD_RK3(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaCDFKD_RK2(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaCDFKD_RK1(cdf, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaCDFKS_RK5(cdf, x, kappa, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaCDFKS_RK4(cdf, x, kappa, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaCDFKS_RK3(cdf, x, kappa, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaCDFKS_RK2(cdf, x, kappa, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaCDFKS_RK1(cdf, x, kappa, invSigma, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaCDFKS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, kappa, invSigma
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar or array of arbitrary rank of Gamma-distributed random values with the specified shape and scale
    !>  parameters \f$(\kappa, \sigma)\f$ of the Gamma distribution corresponding to the procedure arguments `(kappa, sigma)`.
    !>
    !>  See the documentation of [pm_distGamma](@ref pm_distGamma) for more
    !>  information on the Probability Density Function (PDF) of the Gamma distribution.
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
    !>                              On output, it contains Gamma-distributed random value(s).<br>
    !>  \param[in]      kappa   :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `rand`,
    !>                              representing the shape parameter of the Gamma distribution.<br>
    !>  \param[in]      sigma   :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `rand`,
    !>                              representing the scale parameter of the [Gamma distribution](@ref pm_distGamma).<br>
    !>
    !>  \interface{setGammaRand}
    !>  \code{.F90}
    !>
    !>      use pm_distGamma, only: setGammaRand
    !>
    !>      call setGammaRand(rand, kappa, sigma)
    !>      call setGammaRand(rand(..), kappa, sigma)
    !>      call setGammaRand(rng, rand, kappa, sigma)
    !>      call setGammaRand(rng, rand(:), kappa, sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
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
    !>  [getGammaLogPDF](@ref pm_distGamma::getGammaLogPDF)<br>
    !>  [setGammaLogPDF](@ref pm_distGamma::setGammaLogPDF)<br>
    !>  [getGammaCDF](@ref pm_distGamma::getGammaCDF)<br>
    !>  [setGammaCDF](@ref pm_distGamma::setGammaCDF)<br>
    !>
    !>  \example{setGammaRand}
    !>  \include{lineno} example/pm_distGamma/setGammaRand/main.F90
    !>  \compilef{setGammaRand}
    !>  \output{setGammaRand}
    !>  \include{lineno} example/pm_distGamma/setGammaRand/main.out.F90
    !>  \postproc{setGammaRand}
    !>  \include{lineno} example/pm_distGamma/setGammaRand/main.py
    !>  \vis{setGammaRand}
    !>  \image html pm_distGamma/setGammaRand/setGammaRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGamma](@ref test_pm_distGamma)
    !>
    !>  \final{setGammaRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGammaRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGammaRandRNGD_D0_RK5(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGammaRandRNGD_D0_RK4(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGammaRandRNGD_D0_RK3(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGammaRandRNGD_D0_RK2(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGammaRandRNGD_D0_RK1(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGammaRandRNGF_D0_RK5(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGammaRandRNGF_D0_RK4(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGammaRandRNGF_D0_RK3(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGammaRandRNGF_D0_RK2(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGammaRandRNGF_D0_RK1(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGammaRandRNGX_D0_RK5(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGammaRandRNGX_D0_RK4(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGammaRandRNGX_D0_RK3(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGammaRandRNGX_D0_RK2(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGammaRandRNGX_D0_RK1(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setGammaRandRNGD_D1_RK5(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setGammaRandRNGD_D1_RK4(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setGammaRandRNGD_D1_RK3(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setGammaRandRNGD_D1_RK2(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setGammaRandRNGD_D1_RK1(rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setGammaRandRNGF_D1_RK5(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setGammaRandRNGF_D1_RK4(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setGammaRandRNGF_D1_RK3(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setGammaRandRNGF_D1_RK2(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setGammaRandRNGF_D1_RK1(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGF_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGammaRandRNGX_D1_RK5(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGammaRandRNGX_D1_RK4(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGammaRandRNGX_D1_RK3(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGammaRandRNGX_D1_RK2(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGammaRandRNGX_D1_RK1(rng, rand, kappa, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaRandRNGX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: kappa, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distGamma
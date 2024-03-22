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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Raised Cosine distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Raised Cosine distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function **(ICDF)** or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** \f$\pi(\cdot)\f$ of the <b>Raised Cosine distribution</b> over a strictly-positive support \f$x \in [\mu - \sigma, \mu + \sigma]\f$
!>  is defined with the two <b>(location, scale)</b> parameters \f$(\mu, \sigma)\f$ as,<br>
!>  \f{equation}{
!>      \large
!>      \pi(x | \mu, \sigma) = \frac{1}{2\sigma} \left[ 1 + \cos \left({\frac {x - \mu}{\sigma}} ~ \pi \right) \right] ~,~ x \in [\mu - \sigma, \mu + \sigma]
!>  \f}
!>  where \f$\pi\f$ on the right hand side expression is the number \f$\ms{Pi}\f$.<br>
!>
!>  The **CDF** \f$\ms{CDF}(\cdot)\f$ of the <b>Raised Cosine distribution</b> over a strictly-positive support \f$x \in [\mu - \sigma, \mu + \sigma]\f$
!>  is defined with the two <b>(location, scale)</b> parameters \f$(\mu, \sigma)\f$ as,<br>
!>  \f{equation}{
!>      \large
!>      \ms{CDF}(x | \mu, \sigma) = \frac {1}{2} \left[ 1 + \frac{x - \mu}{\sigma} + \frac{1}{\pi} \sin\left( \frac{x - \mu}{\sigma} \pi \right) \right] ~,~ x \in [\mu - \sigma, \mu + \sigma]
!>  \f}
!>  where \f$\pi\f$ on the right hand side expression is the number \f$\ms{Pi}\f$.<br>
!>
!>  \note
!>  The **mean**, **median**, and **mode** of the distribution correspond to the location parameter \f$\mu\f$.<br>
!>
!>  \note
!>  The **variance** of the distribution can be computed as,
!>  \f{equation}{
!>      \large
!>      \ms{variance} = \sigma^2 \left( \frac{1}{3} - \frac{2}{\pi^2} \right) ~.
!>  \f}
!>
!>  \note
!>  The distribution has a constant **skewness** of \f$0\f$.<br>
!>
!>  \note
!>  The distribution has a constant **Excess Kurtosis** of,
!>  \f{equation}{
!>      \large
!>      \ms{Ex. Kurtosis} = \frac{6(90 - \pi^{4})}{5(\pi^{2} - 6)^{2}} \approx -0.59376 ~.
!>  \f}
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>
!>  \test
!>  [test_pm_distCosRaised](@ref test_pm_distCosRaised)
!>
!>  \todo
!>  \pmed
!>  Two additional interfaces for computing the quantiles and random values of Raised Cosine Distribution must be added.<br>
!>  The methodology employed for the [Beta distribution](@ pm_distBeta) might be useful here.<br>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distCosRaised

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distCosRaised"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Probability Density Function (PDF) of the Raised Cosine distribution
    !>  for an input `x` within the support of the distribution \f$x \in [\mu - \sigma, \mu + \sigma]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distCosRaised](@ref pm_distCosRaised) for more information on the Raised Cosine distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the PDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `x`,
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = \f$0.\f$)
    !>  \param[in]  sigma       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                              containing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = \f$1.\f$)
    !>
    !>  \return
    !>  `pdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the PDF of the distribution.<br>
    !>
    !>  \interface{getCosRaisedPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distCosRaised, only: getCosRaisedPDF
    !>
    !>      pdf = getCosRaisedPDF(x, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x >= mu - sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x <= mu + sigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setCosRaisedPDF](@ref pm_distCosRaised::setCosRaisedPDF)<br>
    !>
    !>  \example{getCosRaisedPDF}
    !>  \include{lineno} example/pm_distCosRaised/getCosRaisedPDF/main.F90
    !>  \compilef{getCosRaisedPDF}
    !>  \output{getCosRaisedPDF}
    !>  \include{lineno} example/pm_distCosRaised/getCosRaisedPDF/main.out.F90
    !>  \postproc{getCosRaisedPDF}
    !>  \include{lineno} example/pm_distCosRaised/getCosRaisedPDF/main.py
    !>  \vis{getCosRaisedPDF}
    !>  \image html pm_distCosRaised/getCosRaisedPDF/getCosRaisedPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distCosRaised](@ref test_pm_distCosRaised)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getCosRaisedPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getCosRaisedPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getCosRaisedPDF_RK5(x, mu, sigma) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: pdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getCosRaisedPDF_RK4(x, mu, sigma) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: pdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getCosRaisedPDF_RK3(x, mu, sigma) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: pdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getCosRaisedPDF_RK2(x, mu, sigma) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: pdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getCosRaisedPDF_RK1(x, mu, sigma) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: pdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Probability Density Function (PDF) of the Raised Cosine distribution
    !>  for an input `x` within the support of the distribution \f$x \in [\mu - \sigma, \mu + \sigma]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distCosRaised](@ref pm_distCosRaised) for more information on the Raised Cosine distribution.
    !>
    !>  \param[out] pdf         :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the PDF of the distribution.<br>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the PDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `x`,
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = \f$0.\f$. It must be present if `invSigma` is also present.)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                              containing the inverse of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = \f$1.\f$. It can be present only if `mu` is also present.)
    !>
    !>  \interface{setCosRaisedPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distCosRaised, only: setCosRaisedPDF
    !>
    !>      call setCosRaisedPDF(pdf, x)
    !>      call setCosRaisedPDF(pdf, x, mu)
    !>      call setCosRaisedPDF(pdf, x, mu, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x >= mu - 1 / sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x <= mu + 1 / sigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setCosRaisedPDF](@ref pm_distCosRaised::setCosRaisedPDF)<br>
    !>
    !>  \example{setCosRaisedPDF}
    !>  \include{lineno} example/pm_distCosRaised/setCosRaisedPDF/main.F90
    !>  \compilef{setCosRaisedPDF}
    !>  \output{setCosRaisedPDF}
    !>  \include{lineno} example/pm_distCosRaised/setCosRaisedPDF/main.out.F90
    !>  \postproc{setCosRaisedPDF}
    !>  \include{lineno} example/pm_distCosRaised/setCosRaisedPDF/main.py
    !>  \vis{setCosRaisedPDF}
    !>  \image html pm_distCosRaised/setCosRaisedPDF/setCosRaisedPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distCosRaised](@ref test_pm_distCosRaised)
    !>
    !>  \todo
    !>  \pvlow
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setCosRaisedPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setCosRaisedPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXDD_RK5(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXDD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXDD_RK4(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXDD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXDD_RK3(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXDD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXDD_RK2(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXDD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXDD_RK1(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXDD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMD_RK5(pdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMD_RK4(pdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMD_RK3(pdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMD_RK2(pdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMD_RK1(pdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMI_RK5(pdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMI_RK4(pdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMI_RK3(pdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMI_RK2(pdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setCosRaisedPDFXMI_RK1(pdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedPDFXMI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the Raised Cosine distribution
    !>  for an input `x` within the support of the distribution \f$x \in [\mu - \sigma, \mu + \sigma]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distCosRaised](@ref pm_distCosRaised) for more information on the Raised Cosine distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the CDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `x`,
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = \f$0.\f$)
    !>  \param[in]  sigma       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                              containing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = \f$1.\f$)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution.<br>
    !>
    !>  \interface{getCosRaisedCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distCosRaised, only: getCosRaisedCDF
    !>
    !>      cdf = getCosRaisedCDF(x, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x >= mu - sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x <= mu + sigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setCosRaisedCDF](@ref pm_distCosRaised::setCosRaisedCDF)<br>
    !>
    !>  \example{getCosRaisedCDF}
    !>  \include{lineno} example/pm_distCosRaised/getCosRaisedCDF/main.F90
    !>  \compilef{getCosRaisedCDF}
    !>  \output{getCosRaisedCDF}
    !>  \include{lineno} example/pm_distCosRaised/getCosRaisedCDF/main.out.F90
    !>  \postproc{getCosRaisedCDF}
    !>  \include{lineno} example/pm_distCosRaised/getCosRaisedCDF/main.py
    !>  \vis{getCosRaisedCDF}
    !>  \image html pm_distCosRaised/getCosRaisedCDF/getCosRaisedCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distCosRaised](@ref test_pm_distCosRaised)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getCosRaisedCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getCosRaisedCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getCosRaisedCDF_RK5(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getCosRaisedCDF_RK4(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getCosRaisedCDF_RK3(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getCosRaisedCDF_RK2(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getCosRaisedCDF_RK1(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCosRaisedCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, sigma
        real(RKC)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the Raised Cosine distribution
    !>  for an input `x` within the support of the distribution \f$x \in [\mu - \sigma, \mu + \sigma]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distCosRaised](@ref pm_distCosRaised) for more information on the Raised Cosine distribution.
    !>
    !>  \param[out] cdf         :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution.<br>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the CDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `x`,
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = \f$0.\f$. It must be present if `invSigma` is also present.)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                              containing the inverse of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = \f$1.\f$. It can be present only if `mu` is also present.)
    !>
    !>  \interface{setCosRaisedCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distCosRaised, only: setCosRaisedCDF
    !>
    !>      call setCosRaisedCDF(cdf, x)
    !>      call setCosRaisedCDF(cdf, x, mu)
    !>      call setCosRaisedCDF(cdf, x, mu, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x >= mu - 1 / sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `x <= mu + 1 / sigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setCosRaisedCDF](@ref pm_distCosRaised::setCosRaisedCDF)<br>
    !>
    !>  \example{setCosRaisedCDF}
    !>  \include{lineno} example/pm_distCosRaised/setCosRaisedCDF/main.F90
    !>  \compilef{setCosRaisedCDF}
    !>  \output{setCosRaisedCDF}
    !>  \include{lineno} example/pm_distCosRaised/setCosRaisedCDF/main.out.F90
    !>  \postproc{setCosRaisedCDF}
    !>  \include{lineno} example/pm_distCosRaised/setCosRaisedCDF/main.py
    !>  \vis{setCosRaisedCDF}
    !>  \image html pm_distCosRaised/setCosRaisedCDF/setCosRaisedCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distCosRaised](@ref test_pm_distCosRaised)
    !>
    !>  \todo
    !>  \pvlow
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setCosRaisedCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setCosRaisedCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXDD_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXDD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXDD_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXDD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXDD_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXDD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXDD_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXDD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXDD_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXDD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMD_RK5(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMD_RK4(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMD_RK3(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMD_RK2(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMD_RK1(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMI_RK5(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMI_RK4(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMI_RK3(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMI_RK2(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setCosRaisedCDFXMI_RK1(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCosRaisedCDFXMI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distCosRaised
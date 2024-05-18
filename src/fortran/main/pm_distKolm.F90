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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Kolmogorov distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Kolmogorov distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The Kolmogorov distribution is the distribution of the random variable
!>  \f{equation}{
!>      \large
!>      X = \sup_{t\in [0,1]}|B(t)| ~,
!>  \f}
!>  where \f$B(t)\f$ is the Brownian bridge.<br>
!>  The cumulative distribution function (CDF) of \f$K\f$ over the non-negative support \f$X \in [0, +\infty)\f$ is given by,
!>  \f{equation}{
!>      \large
!>      \ms{CDF}(X\leq x) = 1 - 2\sum_{k=1}^{\infty}(-1)^{k-1}e^{-2k^{2}x^{2}} = {\frac{\sqrt{2\pi}}{x}}\sum_{k=1}^{\infty}e^{-(2k-1)^{2}\pi ^{2}/(8x^{2})} ~,
!>  \f}
!>  which can also be expressed by the Jacobi theta function \f$\vartheta_{01}(z=0;\tau =2ix^{2}/\pi)\f$.<br>
!>  The distribution is named after [Andrey Kolmogorov](https://en.wikipedia.org/wiki/Andrey_Kolmogorov).<br>
!>
!>  The corresponding **PDF** \f$\pi(\cdot)\f$ of the <b>Kolmogorov distribution</b> can be obtained by taking the derivative of the CDF of the distribution with respect to \f$x\f$,
!>  \f{eqnarray}{
!>      \large
!>      \pi(x)
!>      &=& 8x \sum_{k=1}^{\infty}(-1)^{k-1} k^2 e^{-2k^{2}x^{2}} \nonumber \\
!>      &=& {\frac{\sqrt{2\pi}}{x^4}}\sum_{k=1}^{\infty} \left[(2k - 1)^2\pi^2/4 - x^2\right] e^{-(2k-1)^{2}\pi ^{2}/(8x^{2})} \nonumber \\
!>      &=& {\frac{2\sqrt{2\pi}}{x^2}}\sum_{k=1}^{\infty} \left[\frac{(2k - 1)^2\pi^2}{8x^2} - \frac{1}{2}\right] e^{-(2k-1)^{2}\pi ^{2}/(8x^{2})}
!>      ~,
!>  \f}
!>  where the symbol \f$\pi\f$ on the righthand side represents the mathematical number \f$\pi = 3.1415\ldots\f$.<br>
!>
!>  Quantile Function
!>  -----------------
!>
!>  There is no close form expression for the inverse CDF of the Kolmogorov distribution.<br>
!>  However, [root finding methods](@ref pm_mathRoot) can be used to refine an initial guess toward an acceptable answer.<br>
!>
!>  Random Number Generation
!>  ------------------------
!>
!>  In the most naive scenario, the quantile function can be used for random number generation.<br>
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distanceKolm](@ref pm_distanceKolm)<br>
!>
!>  \test
!>  [test_pm_distKolm](@ref test_pm_distKolm)
!>
!>  \todo
!>  \pmed
!>  Two additional interfaces for computing the quantiles and random values of Kolmogorov Distribution must be added.<br>
!>  The methodology employed for the [Beta distribution](@ pm_distBeta) might be useful here.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distKolm

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distKolm"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Kolmogorov
    !>  as defined in the description of [pm_distKolm](@ref pm_distKolm).
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for the definition of the Kolmogorov distribution.
    !>
    !>  \interface{distKolm_type}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: distKolm_type
    !>      type(distKolm_type) :: distKolm
    !>
    !>      distKolm = distKolm_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distKolm_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distKolm_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Probability Density Function (PDF) of the Kolmogorov distribution
    !>  for an input `x` within the support of the distribution \f$X \in [0, +\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the Kolmogorov distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the PDF must be computed.<br>
    !>
    !>  \return
    !>  `pdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the PDF of the distribution.<br>
    !>
    !>  \interface{getKolmPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: getKolmPDF
    !>
    !>      pdf = getKolmPDF(x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. <= x` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setKolmPDF](@ref pm_distKolm::setKolmPDF)<br>
    !>
    !>  \example{getKolmPDF}
    !>  \include{lineno} example/pm_distKolm/getKolmPDF/main.F90
    !>  \compilef{getKolmPDF}
    !>  \output{getKolmPDF}
    !>  \include{lineno} example/pm_distKolm/getKolmPDF/main.out.F90
    !>  \postproc{getKolmPDF}
    !>  \include{lineno} example/pm_distKolm/getKolmPDF/main.py
    !>  \vis{getKolmPDF}
    !>  \image html pm_distKolm/getKolmPDF/getKolmPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getKolmPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getKolmPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getKolmPDF_RK5(x) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: pdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getKolmPDF_RK4(x) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: pdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getKolmPDF_RK3(x) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: pdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getKolmPDF_RK2(x) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: pdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getKolmPDF_RK1(x) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: pdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Probability Density Function (PDF) of the Kolmogorov distribution
    !>  for an input `x` within the support of the distribution \f$X \in [0, +\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the Kolmogorov distribution.
    !>
    !>  \param[out] pdf         :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the PDF of the distribution.<br>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the PDF must be computed.<br>
    !>
    !>  \interface{setKolmPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: setKolmPDF
    !>
    !>      call setKolmPDF(pdf, x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. <= x` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setKolmPDF](@ref pm_distKolm::setKolmPDF)<br>
    !>
    !>  \example{setKolmPDF}
    !>  \include{lineno} example/pm_distKolm/setKolmPDF/main.F90
    !>  \compilef{setKolmPDF}
    !>  \output{setKolmPDF}
    !>  \include{lineno} example/pm_distKolm/setKolmPDF/main.out.F90
    !>  \postproc{setKolmPDF}
    !>  \include{lineno} example/pm_distKolm/setKolmPDF/main.py
    !>  \vis{setKolmPDF}
    !>  \image html pm_distKolm/setKolmPDF/setKolmPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvlow
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setKolmPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setKolmPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setKolmPDF_RK5(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setKolmPDF_RK4(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setKolmPDF_RK3(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setKolmPDF_RK2(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setKolmPDF_RK1(pdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the Kolmogorov distribution
    !>  for an input `x` within the support of the distribution \f$X \in [0, +\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the Kolmogorov distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the CDF must be computed.<br>
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution.<br>
    !>
    !>  \interface{getKolmCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: getKolmCDF
    !>
    !>      cdf = getKolmCDF(x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. <= x` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setKolmCDF](@ref pm_distKolm::setKolmCDF)<br>
    !>
    !>  \example{getKolmCDF}
    !>  \include{lineno} example/pm_distKolm/getKolmCDF/main.F90
    !>  \compilef{getKolmCDF}
    !>  \output{getKolmCDF}
    !>  \include{lineno} example/pm_distKolm/getKolmCDF/main.out.F90
    !>  \postproc{getKolmCDF}
    !>  \include{lineno} example/pm_distKolm/getKolmCDF/main.py
    !>  \vis{getKolmCDF}
    !>  \image html pm_distKolm/getKolmCDF/getKolmCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getKolmCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getKolmCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getKolmCDF_RK5(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getKolmCDF_RK4(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getKolmCDF_RK3(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getKolmCDF_RK2(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getKolmCDF_RK1(x) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x
        real(RKC)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the Kolmogorov distribution
    !>  for an input `x` within the support of the distribution \f$X \in [0, +\infty)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the Kolmogorov distribution.
    !>
    !>  \param[out] cdf         :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution.<br>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments, of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the value at which the CDF must be computed.<br>
    !>
    !>  \interface{setKolmCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: setKolmCDF
    !>
    !>      call setKolmCDF(cdf, x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setKolmCDF](@ref pm_distKolm::setKolmCDF)<br>
    !>
    !>  \example{setKolmCDF}
    !>  \include{lineno} example/pm_distKolm/setKolmCDF/main.F90
    !>  \compilef{setKolmCDF}
    !>  \output{setKolmCDF}
    !>  \include{lineno} example/pm_distKolm/setKolmCDF/main.out.F90
    !>  \postproc{setKolmCDF}
    !>  \include{lineno} example/pm_distKolm/setKolmCDF/main.py
    !>  \vis{setKolmCDF}
    !>  \image html pm_distKolm/setKolmCDF/setKolmCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvlow
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setKolmCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setKolmCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setKolmCDF_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setKolmCDF_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setKolmCDF_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setKolmCDF_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setKolmCDF_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of the quantile corresponding to the specified CDF of <b>Kolmogorov distribution</b>.
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the quantile of the Kolmogorov distribution.<br>
    !>
    !>  \param[in]  cdf :   The input scalar (or array of the same rank, shape, and size as other array like arguments) of,
    !>                      <ol>
    !>                          <li>    type `real` of kind \RKALL,<br>
    !>                      </ol>
    !>                      containing the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>
    !>  \return
    !>  `quan`          :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                      of the same type and kind as the input `cdf`, containing the quantile corresponding to the input `cdf`.<br>
    !>                      By definition, the condition `0 <= quan < +Inf` holds.<br>
    !>
    !>  \interface{getKolmQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: getKolmQuan
    !>
    !>      quan = getKolmQuan(cdf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= cdf < 1` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setKolmQuan](@ref pm_distKolm::setKolmQuan)<br>
    !>
    !>  \example{getKolmQuan}
    !>  \include{lineno} example/pm_distKolm/getKolmQuan/main.F90
    !>  \compilef{getKolmQuan}
    !>  \output{getKolmQuan}
    !>  \include{lineno} example/pm_distKolm/getKolmQuan/main.out.F90
    !>  \postproc{getKolmQuan}
    !>  \include{lineno} example/pm_distKolm/getKolmQuan/main.py
    !>  \vis{getKolmQuan}
    !>  \image html pm_distKolm/getKolmQuan/getKolmQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getKolmQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getKolmQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getKolmQuan_RK5(cdf) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmQuan_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: cdf
        real(RKC)                               :: quan
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getKolmQuan_RK4(cdf) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmQuan_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: cdf
        real(RKC)                               :: quan
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getKolmQuan_RK3(cdf) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmQuan_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: cdf
        real(RKC)                               :: quan
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getKolmQuan_RK2(cdf) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmQuan_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: cdf
        real(RKC)                               :: quan
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getKolmQuan_RK1(cdf) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmQuan_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: cdf
        real(RKC)                               :: quan
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the quantile corresponding to the specified CDF of <b>Kolmogorov distribution</b>.
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the quantile of the Kolmogorov distribution.<br>
    !>
    !>  \param[out] quan        :   The output scalar (or array of the same rank, shape, and size as other array like arguments) of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the quantile corresponding to the input `cdf`.<br>
    !>                              By definition, the condition `0 <= quan < +Inf` holds.<br>
    !>  \param[in]  cdf         :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as the output `quan`,
    !>                              containing the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>
    !>  \interface{setKolmQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: setKolmQuan
    !>
    !>      call setKolmQuan(quan, cdf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= cdf < 1` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getKolmQuan](@ref pm_distKolm::getKolmQuan)<br>
    !>
    !>  \example{setKolmQuan}
    !>  \include{lineno} example/pm_distKolm/setKolmQuan/main.F90
    !>  \compilef{setKolmQuan}
    !>  \output{setKolmQuan}
    !>  \include{lineno} example/pm_distKolm/setKolmQuan/main.out.F90
    !>  \postproc{setKolmQuan}
    !>  \include{lineno} example/pm_distKolm/setKolmQuan/main.py
    !>  \vis{setKolmQuan}
    !>  \image html pm_distKolm/setKolmQuan/setKolmQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \final{setKolmQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setKolmQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setKolmQuan_RK5(quan, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmQuan_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: quan
        real(RKC)   , intent(in)                    :: cdf
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setKolmQuan_RK4(quan, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmQuan_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: quan
        real(RKC)   , intent(in)                    :: cdf
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setKolmQuan_RK3(quan, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmQuan_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: quan
        real(RKC)   , intent(in)                    :: cdf
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setKolmQuan_RK2(quan, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmQuan_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: quan
        real(RKC)   , intent(in)                    :: cdf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setKolmQuan_RK1(quan, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmQuan_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: quan
        real(RKC)   , intent(in)                    :: cdf
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of the random value(s) from the <b>Kolmogorov distribution</b>.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more information on the Kolmogorov distribution.<br>
    !>
    !>  \param[in]  unif    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          containing uniformly-distributed random value(s) in the range \f$[0, 1)\f$.<br>
    !>                          Such uniform random values can be readily obtained via the Fortran intrinsic `random_number()` or
    !>                          [getUnifRand](@ref pm_distUnif::getUnifRand) or [setUnifRand](@ref pm_distUnif::setUnifRand).<br>
    !>
    !>  \return
    !>  `rand`              :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                          of the same type and kind as `unif`, containing the random value(s) from the Kolmogorov distribution.<br>
    !>
    !>  \interface{getKolmRand}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: getKolmRand
    !>
    !>      rand = getKolmRand(unif)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= unif < 1` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setKolmRand](@ref pm_distKolm::setKolmRand)<br>
    !>
    !>  \example{getKolmRand}
    !>  \include{lineno} example/pm_distKolm/getKolmRand/main.F90
    !>  \compilef{getKolmRand}
    !>  \output{getKolmRand}
    !>  \include{lineno} example/pm_distKolm/getKolmRand/main.out.F90
    !>  \postproc{getKolmRand}
    !>  \include{lineno} example/pm_distKolm/getKolmRand/main.py
    !>  \vis{getKolmRand}
    !>  \image html pm_distKolm/getKolmRand/getKolmRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getKolmRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getKolmRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getKolmRand_RK5(unif) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmRand_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: unif
        real(RKC)                               :: rand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getKolmRand_RK4(unif) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmRand_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: unif
        real(RKC)                               :: rand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getKolmRand_RK3(unif) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmRand_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: unif
        real(RKC)                               :: rand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getKolmRand_RK2(unif) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmRand_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: unif
        real(RKC)                               :: rand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getKolmRand_RK1(unif) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKolmRand_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: unif
        real(RKC)                               :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the random value(s) from the
    !>  <b>Kolmogorov distribution</b>.
    !>
    !>  See the documentation of [pm_distKolm](@ref pm_distKolm) for more
    !>  information on generating random numbers from the Kolmogorov distribution.
    !>
    !>  \param[out] rand    :   The input/output scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL.<br>
    !>                          </ol>
    !>                          containing the Kolmogorov-distributed random value(s) in the range \f$[0, +\infty)\f$.<br>
    !>  \param[in]  unif    :   The input scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                          of the same type and kind as the output `rand`, containing uniformly-distributed random value(s) in the range \f$[0, 1)\f$.<br>
    !>
    !>  \interface{setKolmRand}
    !>  \code{.F90}
    !>
    !>      use pm_distKolm, only: setKolmRand
    !>
    !>      call setKolmRand(rand, unif)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= unif < 1` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getKolmRand](@ref pm_distKolm::getKolmRand)<br>
    !>
    !>  \example{setKolmRand}
    !>  \include{lineno} example/pm_distKolm/setKolmRand/main.F90
    !>  \compilef{setKolmRand}
    !>  \output{setKolmRand}
    !>  \include{lineno} example/pm_distKolm/setKolmRand/main.out.F90
    !>  \postproc{setKolmRand}
    !>  \include{lineno} example/pm_distKolm/setKolmRand/main.py
    !>  \vis{setKolmRand}
    !>  \image html pm_distKolm/setKolmRand/setKolmRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distKolm](@ref test_pm_distKolm)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This interface can be extended to support vector-like `rand` arguments other than the `elemental` approach.<br>
    !>  Such an extension would be sensible only if the new interface improves the performance against the `elemental` approach.<br>
    !>
    !>  \final{setKolmRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setKolmRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setKolmRand_RK5(rand, unif)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmRand_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: unif
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setKolmRand_RK4(rand, unif)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmRand_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: unif
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setKolmRand_RK3(rand, unif)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmRand_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: unif
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setKolmRand_RK2(rand, unif)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmRand_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: unif
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setKolmRand_RK1(rand, unif)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKolmRand_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: unif
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distKolm
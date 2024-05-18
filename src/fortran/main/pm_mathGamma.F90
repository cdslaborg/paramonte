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
!>  This module contains procedures and generic interfaces for the Lower and Upper Incomplete Gamma functions.
!>
!>  \details
!>  This module provides multiple function and subroutine procedures for computing the Lower and Upper Incomplete Gamma functions.
!>  These routines mostly differ only in terms of performance and usage convenience.
!>      -#  **If performance is important**, use the subroutine interfaces
!>          [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow) and
!>          [setGammaIncUpp](@ref pm_mathGamma::setGammaIncUpp) to compute
!>          the Lower and Upper Incomplete Gamma functions, respectively.
!>      -#  **If ease of use matters** more than performance, use the function interfaces
!>          [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow) and
!>          [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp) to compute
!>          the Lower and Upper Incomplete Gamma functions, respectively.
!>      -#  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries) and
!>          [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)
!>          are not meant to be used directly to compute directly the incomplete
!>          Gamma functions since their performance depends on the values of the
!>          input parameters.<br>
!>          Rather they are low-level implementations to be called only by the
!>          above higher-level interfaces.<br>
!>          See below for the relevant benchmark.
!>
!>  \benchmarks
!>
!>  \benchmark{setGammaIncLowSeries_vs_setGammaIncUppContFrac, The runtime performance of [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries) vs. [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)}
!>  \include{lineno} benchmark/pm_mathGamma/setGammaIncLowSeries_vs_setGammaIncUppContFrac/main.F90
!>  \compilefb{setGammaIncLowSeries_vs_setGammaIncUppContFrac}
!>  \postprocb{setGammaIncLowSeries_vs_setGammaIncUppContFrac}
!>  \include{lineno} benchmark/pm_mathGamma/setGammaIncLowSeries_vs_setGammaIncUppContFrac/main.py
!>  \visb{setGammaIncLowSeries_vs_setGammaIncUppContFrac}
!>  \image html benchmark/pm_mathGamma/setGammaIncLowSeries_vs_setGammaIncUppContFrac/benchmark.setGammaIncLowSeries_vs_setGammaIncUppContFrac.runtime.png width=1000
!>  \image html benchmark/pm_mathGamma/setGammaIncLowSeries_vs_setGammaIncUppContFrac/benchmark.setGammaIncLowSeries_vs_setGammaIncUppContFrac.runtime.ratio.png width=1000
!>  \moralb{setGammaIncLowSeries_vs_setGammaIncUppContFrac}
!>      -#  The procedures under the generic interface [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)
!>          compute the Lower Incomplete Gamma function suing the series representation of the Incomplete Gamma function while
!>          the procedures under the generic interface [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)
!>          use the Continued Fraction representation of the Gamma function.<br>
!>          The Legendre continued fraction representation is known to converge faster at \f$x > \kappa + 1\f$
!>          whereas the series representation converges faster at \f$x < \kappa + 1\f$.<br>
!>          This performance difference dependence on the specific values of \f$x\f$
!>          and \f$\kappa\f$ is well illustrated by this benchmark.<br>
!>
!>  \test
!>  [test_pm_mathGamma](@ref test_pm_mathGamma)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Monday 12:56 pm, August 16, 2021, Dallas TX

module pm_mathGamma

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathGamma"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **regularized** Lower Incomplete Gamma function for the
    !>  specified shape parameter (\f$\kappa\f$) and upper limit of the integral `x`.
    !>
    !>  \details
    !>  The regularized Lower Incomplete Gamma function is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      P(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_0^{x}~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ should hold, with \f$\kappa\f$ representing the shape parameter of the Gamma function
    !>  (or distribution) and \f$x\f$ representing the upper limit in the integral of the Lower Incomplete Gamma function.<br>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Lower Incomplete Gamma function also represents the *complement* of
    !>  the Cumulative Distribution Function (CDF) of the univariate Gamma distribution with
    !>  the specified shape parameter and standardized `x` (with the scale parameter of unity).
    !>
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series or
    !>                                  its Legendre continued fraction representation.<br>
    !>                                  (**optional**, default is set by [getGammaIncLow()](@ref getGammaIncLow))
    !>
    !>  \return
    !>  `gammaIncLow`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and upper limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge**.
    !>
    !>  \interface{getGammaIncLow}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: getGammaIncLow
    !>
    !>      gammaIncLow = getGammaIncLow(x, kappa, tol = tol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The `kappa` and `x` input arguments must be **positive** `real` numbers.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The procedures under this generic interface solely provide a more convenient method of calling the subroutine
    !>  equivalents [setGammaIncLow](@ref setGammaIncLow).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow)<br>
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow)<br>
    !>  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)<br>
    !>  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)<br>
    !>  Numerical Recipes by Press et al. 1992
    !>
    !>  \example{getGammaIncLow}
    !>  \include{lineno} example/pm_mathGamma/getGammaIncLow/main.F90
    !>  \compilef{getGammaIncLow}
    !>  \output{getGammaIncLow}
    !>  \include{lineno} example/pm_mathGamma/getGammaIncLow/main.out.F90
    !>  \postproc{getGammaIncLow}
    !>  \include{lineno} example/pm_mathGamma/getGammaIncLow/main.py
    !>  \vis{getGammaIncLow}
    !>  \image html pm_mathGamma/getGammaIncLow/getGammaIncLow.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>  \final{getGammaIncLow}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncLow

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncLow_RK5(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncLow
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncLow_RK4(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncLow
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncLow_RK3(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncLow
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncLow_RK2(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncLow
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncLow_RK1(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncLow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **regularized** Upper Incomplete Gamma function for the
    !>  specified shape parameter (\f$\kappa\f$) and lower limit of the integral `x`.
    !>
    !>  \details
    !>  The regularized Upper Incomplete Gamma function is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      Q(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_x^{+\infty}~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ should hold, with \f$\kappa\f$ representing the shape parameter of the Gamma function
    !>  (or distribution) and \f$x\f$ representing the lower limit in the integral of the Upper Incomplete Gamma function.<br>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Upper Incomplete Gamma function also represents the *complement* of
    !>  the Cumulative Distribution Function (CDF) of the univariate Gamma distribution with
    !>  the specified shape parameter and standardized `x` (with the scale parameter of unity).
    !>
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL, representing the lower limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`, representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`, representing the relative accuracy in the convergence checking of the Gamma series or its Legendre continued fraction representation.<br>
    !>                                  (**optional**, default is set by [getGammaIncUpp()](@ref getGammaIncUpp))
    !>
    !>  \return
    !>  `gammaIncUpp`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncUpp` is, by definition, always positive.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge.
    !>
    !>  \interface{getGammaIncUpp}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: getGammaIncUpp
    !>
    !>      gammaIncUpp = getGammaIncUpp(x, kappa, tol = tol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The `kappa` and `x` input arguments must be **positive** `real` numbers.<br>
    !>  Furthermore, `tol << 1.` must hold, if present as an input argument.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The procedures under this generic interface solely provide a more convenient method of calling the subroutine
    !>  equivalents [setGammaIncUpp](@ref setGammaIncUpp).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [setGammaIncUpp](@ref pm_mathGamma::setGammaIncUpp)<br>
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow)<br>
    !>  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)<br>
    !>  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)<br>
    !>
    !>  \example{getGammaIncUpp}
    !>  \include{lineno} example/pm_mathGamma/getGammaIncUpp/main.F90
    !>  \compilef{getGammaIncUpp}
    !>  \output{getGammaIncUpp}
    !>  \include{lineno} example/pm_mathGamma/getGammaIncUpp/main.out.F90
    !>  \postproc{getGammaIncUpp}
    !>  \include{lineno} example/pm_mathGamma/getGammaIncUpp/main.py
    !>  \vis{getGammaIncUpp}
    !>  \image html pm_mathGamma/getGammaIncUpp/getGammaIncUpp.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>  \remark
    !>  See Numerical Recipes by Press et al. 1992 for further details of the Incomplete Gamma function.
    !>
    !>  \final{getGammaIncUpp}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncUpp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncUpp_RK5(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncUpp
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncUpp_RK4(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncUpp
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncUpp_RK3(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncUpp
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncUpp_RK2(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncUpp
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncUpp_RK1(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x, kappa
        real(RKC)   , intent(in)    , optional  :: tol
        real(RKC)                               :: gammaIncUpp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **regularized** Lower Incomplete Gamma function for the
    !>  specified shape parameter (\f$\kappa\f$) and upper limit of the integral `x`.
    !>
    !>  \details
    !>  The regularized Lower Incomplete Gamma function is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      P(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_0^x~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ should hold, with \f$\kappa\f$ representing the shape parameter of the Gamma function
    !>  (or distribution) and \f$x\f$ representing the upper limit in the integral of the Lower Incomplete Gamma function.<br>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Lower Incomplete Gamma function also represents the Cumulative Distribution Function (CDF)
    !>  of the univariate Gamma distribution with the specified shape parameter and standardized `x` (with the scale parameter of unity).
    !>
    !>  \param[out] gammaIncLow     :   The output scalar of the same type and kind as the input argument `x` representing
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and lower limit.
    !>                                  Note that `gammaIncLow` is, by definition, always positive.<br>
    !>  \param[in]  x               :   The input scalar of the type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  logGammaKappa   :   The input scalar or array of the same shape as other input arguments, of the same type and kind
    !>                                  as the input `x`, representing the precomputed \f$\log(\Gamma(\kappa))\f$ which can
    !>                                  be computed by calling the Fortran intrinsic function `log_gamma(kappa)`.<br>
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to (positive) number of iterations taken for the algorithm to converge or its negative if the algorithm fails to converge.<br>
    !>                                  A convergence failure is likely to happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series or its
    !>                                  Legendre continued fraction representation.<br>
    !>                                  (**optional**, if it is missing, the default is set by
    !>                                  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac) or
    !>                                  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)).
    !>
    !>  \interface{setGammaIncLow}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: setGammaIncLow
    !>
    !>      call setGammaIncLow(gammaIncLow, x, logGammaKappa, kappa, info, tol = tol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The `kappa` and `x` input arguments must be **positive** `real` numbers with `logGammaKappa = log_gamma(kappa)` where `log_gamma()` is a Fortran intrinsic function.<br>
    !>  Furthermore, `tol << 1.` must hold, if present as an input argument.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the Gamma and related distributions.<br>
    !>  The logic behind pre-computing and passing `logGammaKappa = log_gamma(kappa)` is to speed up the calculations
    !>  since `log_gamma()` is computationally expensive and its recomputation can be avoided in repeated calls to
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow) with the same shape parameter but different `x`.
    !>
    !>  \see
    !>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp)<br>
    !>  [setGammaIncUpp](@ref pm_mathGamma::setGammaIncUpp)<br>
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)<br>
    !>  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)<br>
    !>  See also The Numerical Recipes by Press et al. 1992 for further details about the Incomplete Gamma function.<br>
    !>
    !>  \example{setGammaIncLow}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncLow/main.F90
    !>  \compilef{setGammaIncLow}
    !>  \output{setGammaIncLow}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncLow/main.out.F90
    !>  \postproc{setGammaIncLow}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncLow/main.py
    !>  \vis
    !>  \image html pm_mathGamma/setGammaIncLow/setGammaIncLow.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>
    !>  \final{setGammaIncLow}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncLow

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncLowDef_RK5(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowDef_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncLowDef_RK4(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowDef_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncLowDef_RK3(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowDef_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncLowDef_RK2(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowDef_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncLowDef_RK1(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowDef_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **regularized** Upper Incomplete Gamma function for the
    !>  specified shape parameter (\f$\kappa\f$) and lower limit of the integral `x`.
    !>
    !>  \details
    !>  The regularized Upper Incomplete Gamma function is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      Q(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_x^{+\infty}~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ should hold, with \f$\kappa\f$ representing the shape parameter of the Gamma function
    !>  (or distribution) and \f$x\f$ representing the lower limit in the integral of the Upper Incomplete Gamma function.<br>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Upper Incomplete Gamma function also represents the *complement*
    !>  of the Cumulative Distribution Function (CDF) of the univariate Gamma distribution
    !>  with the specified shape parameter and standardized `x` (with the scale parameter of unity).
    !>
    !>  \param[out] gammaIncUpp     :   The output scalar of same type and kind as the input argument `x` representing
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and lower limit.
    !>                                  Note that `gammaIncUpp` is, by definition, always positive.<br>
    !>  \param[in]  x               :   The input scalar of the type `real` of kind \RKALL,
    !>                                  representing the lower limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  logGammaKappa   :   The input scalar of the same type and kind
    !>                                  as the input `x`, representing the precomputed \f$\log(\Gamma(\kappa))\f$ which can
    !>                                  be computed by calling the Fortran intrinsic function `log_gamma(kappa)`.<br>
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to (positive) number of iterations taken for the algorithm to converge or its negative if the algorithm fails to converge.<br>
    !>                                  A convergence failure is likely to happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series or its
    !>                                  Legendre continued fraction representation.<br>
    !>                                  (**optional**, if it is missing, the default is set by
    !>                                  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac) or
    !>                                  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)).
    !>
    !>  \interface{setGammaIncUpp}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: setGammaIncUpp
    !>
    !>      call setGammaIncUpp(gammaIncUpp, x, logGammaKappa, kappa, info, tol = tol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The `kappa` and `x` input arguments must be **positive** `real` numbers with `logGammaKappa = log_gamma(kappa)`
    !>  where `log_gamma()` is a Fortran intrinsic function.<br>
    !>  Furthermore, `tol << 1.` must hold, if present as an input argument.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the Gamma and related distributions.<br>
    !>  The logic behind pre-computing and passing `logGammaKappa = log_gamma(kappa)` is to speed up the calculations
    !>  since `log_gamma()` is computationally expensive and its recomputation can be avoided in repeated calls to
    !>  [setGammaIncUpp](@ref pm_mathGamma::setGammaIncUpp) with the same shape parameter but different `x`.
    !>
    !>  \see
    !>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp)<br>
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow)<br>
    !>  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)<br>
    !>  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)<br>
    !>
    !>  \example{setGammaIncUpp}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncUpp/main.F90
    !>  \compilef{setGammaIncUpp}
    !>  \output{setGammaIncUpp}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncUpp/main.out.F90
    !>  \postproc{setGammaIncUpp}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncUpp/main.py
    !>  \vis{setGammaIncUpp}
    !>  \image html pm_mathGamma/setGammaIncUpp/setGammaIncUpp.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>  \remark
    !>  See Numerical Recipes by Press et al. 1992 for further details of the Incomplete Gamma function.
    !>
    !>  \final{setGammaIncUpp}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncUpp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncUppDef_RK5(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppDef_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncUppDef_RK4(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppDef_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncUppDef_RK3(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppDef_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncUppDef_RK2(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppDef_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncUppDef_RK1(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppDef_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **regularized** Lower Incomplete Gamma function for the specified upper limit `x` and shape parameter,
    !>  evaluated by the series representation of the Incomplete Gamma function.
    !>
    !>  \details
    !>  The regularized Lower Incomplete Gamma function is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      P(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_0^{x}~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ are respectively the shape parameter of the Gamma function
    !>  (or distribution) and the upper limit in the integral of the Lower Incomplete Gamma function.<br>
    !>  By definition, the Lower Incomplete Gamma function is always positive.<br>
    !>
    !>  \param[out] gammaIncLow     :   The input scalar of the same type and kind as the input `x`,
    !>                                  representing the regularized Lower Incomplete Gamma function.
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Lower Incomplete Gamma function.
    !>  \param[in]  logGammaKappa   :   The input scalar of the same type and kind as the input `x`,
    !>                                  representing the precomputed \f$\log(\Gamma(\kappa))\f$ which can be computed
    !>                                  by calling the Fortran intrinsic function `log_gamma(kappa)`.<br>
    !>  \param[in]  kappa           :   The input scalar of type `real` of kind \RKALL,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to (positive) number of iterations taken for the series representation to converge or its negative if the series representation fails to converge.<br>
    !>                                  A convergence failure is likely to happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series representation.<br>
    !>                                  (**optional**, default = `10 * epsilon(x)`).
    !>
    !>  \interface{setGammaIncLowSeries}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: setGammaIncLowSeries
    !>
    !>      call setGammaIncLowSeries(gammaIncLow, x, logGammaKappa, kappa, info, tol = tol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The `kappa` and `x` input arguments must be **positive** `real` numbers with `logGammaKappa = log_gamma(kappa)`
    !>  where `log_gamma()` is a Fortran intrinsic function.<br>
    !>  Furthermore, `tol << 1.` must hold, if it is present as an input argument.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the Gamma and related distributions.<br>
    !>  The logic behind pre-computing and passing `logGammaKappa = log_gamma(kappa)` is to speed up the calculations
    !>  since `log_gamma()` is computationally expensive and its recomputation can be avoided in repeated calls to
    !>  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries) with the same shape parameter but different `x` values.
    !>
    !>  \see
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow)<br>
    !>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp)<br>
    !>  [setGammaIncUpp](@ref pm_mathGamma::setGammaIncUpp)<br>
    !>  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac)<br>
    !>  See also The Numerical Recipes by Press et al. 1992 for further details about the Incomplete Gamma function.<br>
    !>
    !>  \example{setGammaIncLowSeries}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncLowSeries/main.F90
    !>  \compilef{setGammaIncLowSeries}
    !>  \output{setGammaIncLowSeries}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncLowSeries/main.out.F90
    !>  \postproc{setGammaIncLowSeries}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncLowSeries/main.py
    !>  \vis{setGammaIncLowSeries}
    !>  \image html pm_mathGamma/setGammaIncLowSeries/setGammaIncLowSeries.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>  \final{setGammaIncLowSeries}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncLowSeries

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncLowSeries_RK5(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeries_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncLowSeries_RK4(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeries_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncLowSeries_RK3(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeries_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncLowSeries_RK2(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeries_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncLowSeries_RK1(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeries_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)               :: gammaIncLow
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **regularized** Upper Incomplete Gamma function for the specified lower limit `x` and shape parameter,
    !>  evaluated by the Legendre continued fraction representation of the Incomplete Gamma function.
    !>
    !>  \details
    !>  The regularized Upper Incomplete Gamma function is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      Q(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_x^{+\infty}~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ are respectively the shape parameter of the Gamma function
    !>  (or distribution) and the lower limit in the integral of the Upper Incomplete Gamma function.<br>
    !>  By definition, the Upper Incomplete Gamma function is always positive.<br>
    !>
    !>  \param[out] gammaIncUpp     :   The input scalar of the same type and kind as the input `x`,
    !>                                  representing the regularized Upper Incomplete Gamma function.
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL,
    !>                                  representing the lower limit in the integral of the Upper Incomplete Gamma function.
    !>  \param[in]  logGammaKappa   :   The input scalar the same type and kind as the input `x`,
    !>                                  representing the precomputed \f$\log(\Gamma(\kappa))\f$ which can be computed
    !>                                  by calling the Fortran intrinsic function `log_gamma(kappa)`.<br>
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as the input `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Upper Incomplete Gamma function.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to (positive) number of iterations taken for the Legendre continued fraction representation to converge
    !>                                  or its negative if the Legendre continued fraction representation fails to converge.<br>
    !>                                  A convergence failure is likely to happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Legendre continued fraction
    !>                                  representation of the Gamma function.<br>
    !>                                  (**optional**, default = `10 * epsilon(x)`).
    !>
    !>  \interface{setGammaIncUppContFrac}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: setGammaIncUppContFrac
    !>
    !>      gammaIncLow = setGammaIncUppContFrac(gammaIncUpp, x, logGammaKappa, kappa, info, tol = tol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The `kappa` and `x` input arguments must be **positive** `real` numbers
    !>  with `logGammaKappa = log_gamma(kappa)` where `log_gamma()` is a Fortran intrinsic function.<br>
    !>  Furthermore, `tol << 1.` must hold, if it is present as an input argument.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the Gamma and related distributions.<br>
    !>  The logic behind pre-computing and passing `logGammaKappa = log_gamma(kappa)` is to speed up the calculations
    !>  since `log_gamma()` is computationally expensive and its recomputation can be avoided in repeated calls to
    !>  [setGammaIncUppContFrac](@ref pm_mathGamma::setGammaIncUppContFrac) with the same shape parameter but different `x` values.
    !>
    !>  \see
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaIncLow](@ref pm_mathGamma::setGammaIncLow)<br>
    !>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp)<br>
    !>  [setGammaIncUpp](@ref pm_mathGamma::setGammaIncUpp)<br>
    !>  [setGammaIncLowSeries](@ref pm_mathGamma::setGammaIncLowSeries)<br>
    !>  See also The Numerical Recipes by Press et al. 1992 for further details about the Incomplete Gamma function.<br>
    !>
    !>  \example{setGammaIncUppContFrac}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncUppContFrac/main.F90
    !>  \compilef{setGammaIncUppContFrac}
    !>  \output{setGammaIncUppContFrac}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncUppContFrac/main.out.F90
    !>  \postproc{setGammaIncUppContFrac}
    !>  \include{lineno} example/pm_mathGamma/setGammaIncUppContFrac/main.py
    !>  \vis{setGammaIncUppContFrac}
    !>  \image html pm_mathGamma/setGammaIncUppContFrac/setGammaIncUppContFrac.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>  \final{setGammaIncUppContFrac}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncUppContFrac

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncUppContFrac_RK5(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFrac_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncUppContFrac_RK4(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFrac_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncUppContFrac_RK3(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFrac_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncUppContFrac_RK2(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFrac_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncUppContFrac_RK1(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFrac_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)               :: gammaIncUpp
        real(RKC)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKC)   , intent(in)    , optional  :: tol
    end subroutine
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !> Return the Gamma function for a half-integer input as real of kind \RK.
!    !>
!    !> \param[in]   positiveHalfInteger :   The input half integer as a real number
!    !>
!    !> \return
!    !> `gammaHalfInt`                   :   The Gamma function for a half integer input.
!    !>
!    !> \remark
!    !> The equation for half-integer Gamma-function is given as,
!    !> \f{equation}{
!    !>     \Gamma \left( \frac{n}{2} \right) = \sqrt \pi \frac{ (n-2)!! }{ 2^\frac{n-1}{2} } ~,
!    !> \f}
!    pure function getGamHalfInt(positiveHalfInteger) result(gammaHalfInt)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGamHalfInt
!#endif
!        use pm_kind, only: IK, RK, SQRT_PI
!        implicit none
!        real(RK), intent(in) :: positiveHalfInteger
!        real(RKC)            :: gammaHalfInt
!        integer(IK)          :: i,k
!        gammaHalfInt = SQRT_PI
!        k = nint(positiveHalfInteger-0.5_RK,kind=IK) ! positiveHalfInteger = k + 1/2
!        do i = k+1, 2*k
!            gammaHalfInt = gammaHalfInt * 0.25_RK * i
!        end do
!    end function getGamHalfInt
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> Return the natural logarithm of the Gamma function for a half-integer input as real of kind \RK.
!    !>
!    !> \param[in]   positiveHalfInteger :   The input half integer as a real number
!    !>
!    !> \return
!    !> `gammaHalfInt`                   :   The Gamma function for a half integer input.
!    !>
!    !> \remark
!    !> The equation for half-integer Gamma-function is given as,
!    !> \f{equation}{
!    !>     \Gamma \left( \frac{n}{2} \right) = \sqrt \pi \frac{ (n-2)!! }{ 2^\frac{n-1}{2} } ~,
!    !> \f}
!    pure function getLogGammaHalfInt(positiveHalfInteger) result(logGammaHalfInt)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogGammaHalfInt
!#endif
!        use pm_kind, only: IK, RK, SQRT_PI
!        implicit none
!        real(RK), intent(in)    :: positiveHalfInteger
!        real(RK), parameter     :: COEF = log(0.25_RK)
!        real(RK), parameter     :: LOG_SQRTPI = log(SQRT_PI)
!        real(RKC)               :: logGammaHalfInt
!        integer(IK)             :: i, k
!        k = nint(positiveHalfInteger-0.5_RK,kind=IK) ! positiveHalfInteger = k + 1/2
!        logGammaHalfInt = LOG_SQRTPI
!        do i = k+1, 2*k
!            logGammaHalfInt = logGammaHalfInt + COEF + log(real(i,kind=RK))
!        end do
!    end function getLogGammaHalfInt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathGamma ! LCOV_EXCL_LINE
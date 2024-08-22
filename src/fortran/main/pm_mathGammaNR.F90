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
!>  These routines mostly differ only in terms of performance and usage convenience.<br>
!>      -#  **If performance is important**, use the subroutine interface
!>          [setGammaIncLowNR](@ref pm_mathGammaNR::setGammaIncLowNR) and
!>          [setGammaIncUppNR](@ref pm_mathGammaNR::setGammaIncUppNR)
!>          to compute the Lower and Upper Incomplete Gamma function.<br>
!>      -#  **If ease of use matters** more than performance, use the function interfaces
!>          [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR) and
!>          [getGammaIncUppNR](@ref pm_mathGammaNR::getGammaIncUppNR)
!>          to compute the Lower and Upper Incomplete Gamma function.<br>
!>      -#  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR) and
!>          [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR)
!>          are not meant to be used directly to compute directly the incomplete
!>          Gamma functions since their performance depends on the values of the input parameters.<br>
!>          Rather they are low-level implementations to be called only by the above higher-level interfaces.<br>
!>          See below for the relevant benchmark.<br>
!>
!>  \warning
!>  Although all generic interfaces of this module are available for  all processor `real` kinds,
!>  the accuracy and performance of the implemented algorithms are optimized for IEEE double precision.<br>
!>  In particular, the algorithms **may not accurately** compute the lower incomplete gamma function
!>  in extended precision (e.g., 128 bits) mode corresponding to \RKH kind type parameter.<br>
!>
!>  \note
!>  The computations of this module are an implementation of the algorithm proposed approach by:<br>
!>  Numerical Recipes in Fortran 77: The Art of Scientific Computing 2nd Edition
!>  William H. Press (Author), Brian P. Flannery (Author), Saul A. Teukolsky (Author), William T. Vetterling (Author)
!>
!>  \benchmarks
!>
!>  \benchmark{setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR, The runtime performance of [setGammaIncLowSeries](@ref pm_mathGammaNR::setGammaIncLowSeriesNR) vs. [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR)}
!>  \include{lineno} benchmark/pm_mathGammaNR/setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR/main.F90
!>  \compilefb{setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR}
!>  \postprocb{setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR}
!>  \include{lineno} benchmark/pm_mathGammaNR/setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR/main.py
!>  \visb{setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR}
!>  \image html benchmark/pm_mathGammaNR/setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR/benchmark.setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR.runtime.png width=1000
!>  \image html benchmark/pm_mathGammaNR/setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR/benchmark.setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR.runtime.ratio.png width=1000
!>  \moralb{setGammaIncLowSeriesNR_vs_setGammaIncUppContFracNR}
!>      -#  The procedures under the generic interface [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR)
!>          compute the Lower Incomplete Gamma function suing the series representation of the Incomplete Gamma function while
!>          the procedures under the generic interface [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR)
!>          use the Continued Fraction representation of the Gamma function.<br>
!>          The Legendre continued fraction representation is known to converge faster at \f$x > \kappa + 1\f$
!>          whereas the series representation converges faster at \f$x < \kappa + 1\f$.<br>
!>          This performance difference dependence on the specific values of \f$x\f$
!>          and \f$\kappa\f$ is well illustrated by this benchmark.<br>
!>
!>  \see
!>  [pm_mathGammaNR](@ref pm_mathGammaNR) for detailed description of the (Regularized Incomplete) Gamma Function.<br>
!>
!>  \test
!>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
!>
!>  \todo
!>  \pvhigh
!>  The implementation of the algorithms of this module must be properly changed to
!>  allow reliable extended-precision computations of the incomplete Gamma function.<br>
!>  This would require significant investment in making the original algorithms of Gil et al. kind-agnostic.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, July 22, 2024, 11:45 AM, NASA Goddard Space Flight Center, Washington, D.C.<br>

module pm_mathGammaNR

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathGammaNR"

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
    !>                                  (**optional**, default is set by [setGammaIncLowNR()](@ref pm_mathGammaNR::setGammaIncLowNR))
    !>
    !>  \return
    !>  `gammaIncLow`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and upper limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge**.
    !>
    !>  \interface{getGammaIncLowNR}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaNR, only: getGammaIncLowNR
    !>
    !>      gammaIncLow = getGammaIncLowNR(x, kappa, tol = tol)
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
    !>  equivalents [setGammaIncLowNR](@ref setGammaIncLowNR).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [setGammaIncLowNR](@ref pm_mathGammaNR::setGammaIncLowNR)<br>
    !>  [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR)<br>
    !>  [setGammaIncLowNR](@ref pm_mathGammaNR::setGammaIncLowNR)<br>
    !>  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR)<br>
    !>  [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR)<br>
    !>  Numerical Recipes by Press et al. 1992
    !>
    !>  \example{getGammaIncLowNR}
    !>  \include{lineno} example/pm_mathGammaNR/getGammaIncLowNR/main.F90
    !>  \compilef{getGammaIncLowNR}
    !>  \output{getGammaIncLowNR}
    !>  \include{lineno} example/pm_mathGammaNR/getGammaIncLowNR/main.out.F90
    !>  \postproc{getGammaIncLowNR}
    !>  \include{lineno} example/pm_mathGammaNR/getGammaIncLowNR/main.py
    !>  \vis{getGammaIncLowNR}
    !>  \image html pm_mathGammaNR/getGammaIncLowNR/getGammaIncLowNR.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
    !>
    !>  \final{getGammaIncLowNR}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncLowNR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncLowNR_RK5(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowNR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncLowNR_RK4(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowNR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncLowNR_RK3(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowNR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncLowNR_RK2(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowNR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncLowNR_RK1(x, kappa, tol) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowNR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncLow
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
    !>                                  (**optional**, default is set by [setGammaIncUppNR()](@ref pm_mathGammaNR::setGammaIncUppNR))
    !>
    !>  \return
    !>  `gammaIncUpp`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncUpp` is, by definition, always positive.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge.
    !>
    !>  \interface{getGammaIncUppNR}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaNR, only: getGammaIncUppNR
    !>
    !>      gammaIncUpp = getGammaIncUppNR(x, kappa, tol = tol)
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
    !>  equivalents [setGammaIncUppNR](@ref setGammaIncUppNR).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR)<br>
    !>  [setGammaIncUppNR](@ref pm_mathGammaNR::setGammaIncUppNR)<br>
    !>
    !>  \example{getGammaIncUppNR}
    !>  \include{lineno} example/pm_mathGammaNR/getGammaIncUppNR/main.F90
    !>  \compilef{getGammaIncUppNR}
    !>  \output{getGammaIncUppNR}
    !>  \include{lineno} example/pm_mathGammaNR/getGammaIncUppNR/main.out.F90
    !>  \postproc{getGammaIncUppNR}
    !>  \include{lineno} example/pm_mathGammaNR/getGammaIncUppNR/main.py
    !>  \vis{getGammaIncUppNR}
    !>  \image html pm_mathGammaNR/getGammaIncUppNR/getGammaIncUppNR.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
    !>
    !>  \remark
    !>  See Numerical Recipes by Press et al. 1992 for further details of the Incomplete Gamma function.
    !>
    !>  \final{getGammaIncUppNR}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncUppNR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncUppNR_RK5(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppNR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncUppNR_RK4(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppNR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncUppNR_RK3(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppNR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncUppNR_RK2(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppNR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncUppNR_RK1(x, kappa, tol) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppNR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)   , intent(in)    , optional  :: tol
        real(RKG)                               :: gammaIncUpp
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
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>  \param[in]  x               :   The input scalar of the type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  logGammaKappa   :   The input scalar or array of the same shape as other input arguments, of the same type and kind
    !>                                  as the input `x`, representing the precomputed \f$\log(\Gamma(\kappa))\f$ which can
    !>                                  be computed by calling the Fortran intrinsic function `log_gamma(kappa)`.<br>
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to (positive) number of iterations taken for the algorithm to converge or its negative if the algorithm fails to converge.<br>
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series or its
    !>                                  Legendre continued fraction representation.<br>
    !>                                  (**optional**, if it is missing, the default is set by
    !>                                  [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR) or
    !>                                  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR)).
    !>
    !>  \interface{setGammaIncLowNR}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaNR, only: setGammaIncLowNR
    !>
    !>      call setGammaIncLowNR(gammaIncLow, x, logGammaKappa, kappa, info, tol = tol)
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
    !>  [setGammaIncLowNR](@ref pm_mathGammaNR::setGammaIncLowNR) with the same shape parameter but different `x`.
    !>
    !>  \see
    !>  [getGammaIncUppNR](@ref pm_mathGammaNR::getGammaIncUppNR)<br>
    !>  [setGammaIncUppNR](@ref pm_mathGammaNR::setGammaIncUppNR)<br>
    !>  [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR)<br>
    !>  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR)<br>
    !>  [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR)<br>
    !>  See also The Numerical Recipes by Press et al. 1992 for further details about the Incomplete Gamma function.<br>
    !>
    !>  \example{setGammaIncLowNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncLowNR/main.F90
    !>  \compilef{setGammaIncLowNR}
    !>  \output{setGammaIncLowNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncLowNR/main.out.F90
    !>  \postproc{setGammaIncLowNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncLowNR/main.py
    !>  \vis
    !>  \image html pm_mathGammaNR/setGammaIncLowNR/setGammaIncLowNR.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
    !>
    !>
    !>  \final{setGammaIncLowNR}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncLowNR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncLowNR_RK5(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowNR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncLowNR_RK4(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowNR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncLowNR_RK3(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowNR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncLowNR_RK2(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowNR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncLowNR_RK1(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowNR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
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
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and lower limit.<br>
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
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series or its
    !>                                  Legendre continued fraction representation.<br>
    !>                                  (**optional**, if it is missing, the default is set by
    !>                                  [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR) or
    !>                                  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR)).
    !>
    !>  \interface{setGammaIncUppNR}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaNR, only: setGammaIncUppNR
    !>
    !>      call setGammaIncUppNR(gammaIncUpp, x, logGammaKappa, kappa, info, tol = tol)
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
    !>  [setGammaIncUppNR](@ref pm_mathGammaNR::setGammaIncUppNR) with the same shape parameter but different `x`.
    !>
    !>  \see
    !>  [getGammaIncUppNR](@ref pm_mathGammaNR::getGammaIncUppNR)<br>
    !>  [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR)<br>
    !>
    !>  \example{setGammaIncUppNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncUppNR/main.F90
    !>  \compilef{setGammaIncUppNR}
    !>  \output{setGammaIncUppNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncUppNR/main.out.F90
    !>  \postproc{setGammaIncUppNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncUppNR/main.py
    !>  \vis{setGammaIncUppNR}
    !>  \image html pm_mathGammaNR/setGammaIncUppNR/setGammaIncUppNR.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
    !>
    !>  \remark
    !>  See Numerical Recipes by Press et al. 1992 for further details of the Incomplete Gamma function.
    !>
    !>  \final{setGammaIncUppNR}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncUppNR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncUppNR_RK5(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppNR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncUppNR_RK4(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppNR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncUppNR_RK3(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppNR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncUppNR_RK2(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppNR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncUppNR_RK1(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppNR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
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
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Gamma series representation.<br>
    !>                                  (**optional**, default = `10 * epsilon(x)`).
    !>
    !>  \interface{setGammaIncLowSeriesNR}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaNR, only: setGammaIncLowSeriesNR
    !>
    !>      call setGammaIncLowSeriesNR(gammaIncLow, x, logGammaKappa, kappa, info, tol = tol)
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
    !>  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR) with the same shape parameter but different `x` values.
    !>
    !>  \see
    !>  [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR)<br>
    !>  [getGammaIncUppNR](@ref pm_mathGammaNR::getGammaIncUppNR)<br>
    !>  [setGammaIncLowNR](@ref pm_mathGammaNR::setGammaIncLowNR)<br>
    !>  [setGammaIncUppNR](@ref pm_mathGammaNR::setGammaIncUppNR)<br>
    !>  [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR)<br>
    !>  See also The Numerical Recipes by Press et al. 1992 for further details about the Incomplete Gamma function.<br>
    !>
    !>  \example{setGammaIncLowSeriesNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncLowSeriesNR/main.F90
    !>  \compilef{setGammaIncLowSeriesNR}
    !>  \output{setGammaIncLowSeriesNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncLowSeriesNR/main.out.F90
    !>  \postproc{setGammaIncLowSeriesNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncLowSeriesNR/main.py
    !>  \vis{setGammaIncLowSeriesNR}
    !>  \image html pm_mathGammaNR/setGammaIncLowSeriesNR/setGammaIncLowSeriesNR.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
    !>
    !>  \final{setGammaIncLowSeriesNR}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncLowSeriesNR

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncLowSeriesNR_RK5(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeriesNR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncLowSeriesNR_RK4(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeriesNR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncLowSeriesNR_RK3(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeriesNR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncLowSeriesNR_RK2(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeriesNR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncLowSeriesNR_RK1(gammaIncLow, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowSeriesNR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
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
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>  \param[in]  tol             :   The input scalar of the same type and kind as `x`,
    !>                                  representing the relative accuracy in the convergence checking of the Legendre continued fraction
    !>                                  representation of the Gamma function.<br>
    !>                                  (**optional**, default = `10 * epsilon(x)`).
    !>
    !>  \interface{setGammaIncUppContFracNR}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaNR, only: setGammaIncUppContFracNR
    !>
    !>      gammaIncLow = setGammaIncUppContFracNR(gammaIncUpp, x, logGammaKappa, kappa, info, tol = tol)
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
    !>  [setGammaIncUppContFracNR](@ref pm_mathGammaNR::setGammaIncUppContFracNR) with the same shape parameter but different `x` values.
    !>
    !>  \see
    !>  [getGammaIncLowNR](@ref pm_mathGammaNR::getGammaIncLowNR)<br>
    !>  [getGammaIncUppNR](@ref pm_mathGammaNR::getGammaIncUppNR)<br>
    !>  [setGammaIncLowNR](@ref pm_mathGammaNR::setGammaIncLowNR)<br>
    !>  [setGammaIncUppNR](@ref pm_mathGammaNR::setGammaIncUppNR)<br>
    !>  [setGammaIncLowSeriesNR](@ref pm_mathGammaNR::setGammaIncLowSeriesNR)<br>
    !>  See also The Numerical Recipes by Press et al. 1992 for further details about the Incomplete Gamma function.<br>
    !>
    !>  \example{setGammaIncUppContFracNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncUppContFracNR/main.F90
    !>  \compilef{setGammaIncUppContFracNR}
    !>  \output{setGammaIncUppContFracNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncUppContFracNR/main.out.F90
    !>  \postproc{setGammaIncUppContFracNR}
    !>  \include{lineno} example/pm_mathGammaNR/setGammaIncUppContFracNR/main.py
    !>  \vis{setGammaIncUppContFracNR}
    !>  \image html pm_mathGammaNR/setGammaIncUppContFracNR/setGammaIncUppContFracNR.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaNR](@ref test_pm_mathGammaNR)
    !>
    !>  \final{setGammaIncUppContFracNR}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncUppContFracNR

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncUppContFracNR_RK5(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFracNR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncUppContFracNR_RK4(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFracNR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncUppContFracNR_RK3(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFracNR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncUppContFracNR_RK2(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFracNR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncUppContFracNR_RK1(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppContFracNR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
        real(RKG)   , intent(in)    , optional  :: tol
    end subroutine
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathGammaNR ! LCOV_EXCL_LINE
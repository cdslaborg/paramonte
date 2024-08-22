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
!>  The current implementation of the functionalities of this module rely on [pm_mathGammaGil](@ref pm_mathGammaGil)
!>  which seems to provide the most reliable estimation of the incomplete Gamma functions for a wide range of input arguments.<br>
!>
!>  \see
!>  [pm_mathGammaAM](@ref pm_mathGammaAM)<br>
!>  [pm_mathGammaNR](@ref pm_mathGammaNR)<br>
!>  [pm_mathGammaGil](@ref pm_mathGammaGil)<br>
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
    !>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>
    !>  The regularized Lower Incomplete Gamma function also represents the *complement* of
    !>  the Cumulative Distribution Function (CDF) of the univariate Gamma distribution with
    !>  the specified shape parameter and standardized `x` (with the scale parameter of unity).<br>
    !>
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>
    !>  \return
    !>  `gammaIncLow`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and upper limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge**.
    !>
    !>  \interface{getGammaIncLow}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: getGammaIncLow
    !>
    !>      gammaIncLow = getGammaIncLow(x, kappa)
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
    !>  equivalents [setGammaInc](@ref setGammaInc).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [setGammaInc](@ref pm_mathGamma::setGammaInc)<br>
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [setGammaInc](@ref pm_mathGamma::setGammaInc)<br>
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
    impure elemental module function getGammaIncLow_RK5(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncLow_RK4(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncLow_RK3(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncLow_RK2(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncLow_RK1(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLow_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x, kappa
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
    !>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>
    !>  The regularized Upper Incomplete Gamma function also represents the *complement* of
    !>  the Cumulative Distribution Function (CDF) of the univariate Gamma distribution with
    !>  the specified shape parameter and standardized `x` (with the scale parameter of unity).<br>
    !>
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL, representing the lower limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`, representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
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
    !>      gammaIncUpp = getGammaIncUpp(x, kappa)
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
    !>  equivalents [setGammaInc](@ref pm_mathGamma::setGammaInc).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp)<br>
    !>  [setGammaInc](@ref pm_mathGamma::setGammaInc)<br>
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
    !>  \final{getGammaIncUpp}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncUpp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncUpp_RK5(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncUpp_RK4(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncUpp_RK3(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncUpp_RK2(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncUpp_RK1(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUpp_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x, kappa
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
    !>  The regularized **Lower Incomplete Gamma function** is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      P(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_0^x~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ should hold, with \f$\kappa\f$ representing the shape parameter of the Gamma function
    !>  (or distribution) and \f$x\f$ representing the upper limit in the integral of the Lower Incomplete Gamma function.<br>
    !>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>
    !>  The regularized Lower Incomplete Gamma function also represents the Cumulative Distribution Function (CDF)
    !>  of the univariate Gamma distribution with the specified shape parameter and standardized `x` (with the scale parameter of unity).<br>
    !>
    !>  The regularized **Upper Incomplete Gamma function** is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      Q(\kappa, x) = \frac{1}{\Gamma(\kappa)} \int_x^{+\infty}~t^{\kappa-1}{\mathrm e}^{-t} ~ dt ~,
    !>  \f}
    !>
    !>  where \f$(\kappa > 0, x > 0)\f$ should hold, with \f$\kappa\f$ representing the shape parameter of the Gamma function
    !>  (or distribution) and \f$x\f$ representing the lower limit in the integral of the Upper Incomplete Gamma function.<br>
    !>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>
    !>  The regularized Upper Incomplete Gamma function also represents the *complement*
    !>  of the Cumulative Distribution Function (CDF) of the univariate Gamma distribution
    !>  with the specified shape parameter and standardized `x` (with the scale parameter of unity).<br>
    !>
    !>  \param[out] gammaIncLow     :   The output scalar of the same type and kind as the input argument `x` representing
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>  \param[out] gammaIncUpp     :   The output scalar of same type and kind as the input argument `x` representing
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncUpp` is, by definition, always positive.<br>
    !>  \param[in]  x               :   The input scalar of the type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to zero if the algorithm succeeds to converge or a negative value if the algorithm fails to converge.<br>
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>
    !>  \interface{setGammaInc}
    !>  \code{.F90}
    !>
    !>      use pm_mathGamma, only: setGammaInc
    !>
    !>      call setGammaInc(gammaIncLow, gammaIncUpp, x, kappa, info)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGammaInc](@ref pm_mathGamma::setGammaInc)<br>
    !>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp)<br>
    !>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow)<br>
    !>
    !>  \example{setGammaInc}
    !>  \include{lineno} example/pm_mathGamma/setGammaInc/main.F90
    !>  \compilef{setGammaInc}
    !>  \output{setGammaInc}
    !>  \include{lineno} example/pm_mathGamma/setGammaInc/main.out.F90
    !>  \postproc{setGammaInc}
    !>  \include{lineno} example/pm_mathGamma/setGammaInc/main.py
    !>  \vis
    !>  \image html pm_mathGamma/setGammaInc/setGammaInc.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGamma](@ref test_pm_mathGamma)
    !>
    !>
    !>  \final{setGammaInc}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaInc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaInc_RK5(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaInc_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaInc_RK4(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaInc_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaInc_RK3(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaInc_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaInc_RK2(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaInc_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaInc_RK1(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaInc_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
!        real(RKG)            :: gammaHalfInt
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
!        real(RKG)               :: logGammaHalfInt
!        integer(IK)             :: i, k
!        k = nint(positiveHalfInteger-0.5_RK,kind=IK) ! positiveHalfInteger = k + 1/2
!        logGammaHalfInt = LOG_SQRTPI
!        do i = k+1, 2*k
!            logGammaHalfInt = logGammaHalfInt + COEF + log(real(i,kind=RK))
!        end do
!    end function getLogGammaHalfInt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathGamma ! LCOV_EXCL_LINE
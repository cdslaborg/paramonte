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
!>          [setGammaIncGil](@ref pm_mathGammaGil::setGammaIncGil)
!>          to compute the Lower and Upper Incomplete Gamma functions.<br>
!>      -#  **If ease of use matters** more than performance, use the function
!>          interfaces [getGammaIncLowGil](@ref pm_mathGammaGil::getGammaIncLowGil)
!>          to compute the Lower and Upper Incomplete Gamma functions.<br>
!>
!>  \warning
!>  Although all generic interfaces of this module are available for  all processor `real` kinds,
!>  the accuracy and performance of the implemented algorithms are optimized for IEEE double precision.<br>
!>  In particular, the algorithms **may not accurately** compute the Lower and Upper Incomplete Gamma functions
!>  in extended precision (e.g., 128 bits) mode corresponding to \RKH kind type parameter.<br>
!>
!>  \note
!>  The computations of this module are explicitly based on the proposed approach by:<br>
!>  Gil et al, 2012, EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE GAMMA FUNCTION RATIOS<br>
!>
!>  \see
!>  [pm_mathGamma](@ref pm_mathGamma) for detailed description of the (Regularized Incomplete) Gamma Function.<br>
!>
!>  \test
!>  [test_pm_mathGammaGil](@ref test_pm_mathGammaGil)
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

module pm_mathGammaGil

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathGammaGil"

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
    !>
    !>  \return
    !>  `gammaIncLow`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Lower Incomplete Gamma function for the specified `kappa` and upper limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge**.
    !>
    !>  \interface{getGammaIncLowGil}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaGil, only: getGammaIncLowGil
    !>
    !>      gammaIncLow = getGammaIncLowGil(x, kappa)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGammaIncLowGil](@ref pm_mathGammaGil::getGammaIncLowGil)<br>
    !>  [setGammaIncGil](@ref pm_mathGammaGil::setGammaIncGil)<br>
    !>
    !>  \example{getGammaIncLowGil}
    !>  \include{lineno} example/pm_mathGammaGil/getGammaIncLowGil/main.F90
    !>  \compilef{getGammaIncLowGil}
    !>  \output{getGammaIncLowGil}
    !>  \include{lineno} example/pm_mathGammaGil/getGammaIncLowGil/main.out.F90
    !>  \postproc{getGammaIncLowGil}
    !>  \include{lineno} example/pm_mathGammaGil/getGammaIncLowGil/main.py
    !>  \vis{getGammaIncLowGil}
    !>  \image html pm_mathGammaGil/getGammaIncLowGil/getGammaIncLowGil.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaGil](@ref test_pm_mathGammaGil)
    !>
    !>  \final{getGammaIncLowGil}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncLowGil

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGammaIncLowGil_RK5(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowGil_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGammaIncLowGil_RK4(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowGil_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGammaIncLowGil_RK3(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowGil_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGammaIncLowGil_RK2(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowGil_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGammaIncLowGil_RK1(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowGil_RK1
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
    !>  specified shape parameter (\f$\kappa\f$) and upper limit of the integral `x`.
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
    !>  (or distribution) and \f$x\f$ representing the upper limit in the integral of the Upper Incomplete Gamma function.<br>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Upper Incomplete Gamma function also represents the *complement* of
    !>  the Cumulative Distribution Function (CDF) of the univariate Gamma distribution with
    !>  the specified shape parameter and standardized `x` (with the scale parameter of unity).<br>
    !>
    !>  \param[in]  x               :   The input scalar of type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Upper Incomplete Gamma function \f$Q(\kappa,x)\f$.<br>
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Upper Incomplete Gamma function \f$Q(\kappa,x)\f$.<br>
    !>
    !>  \return
    !>  `gammaIncUpp`               :   The output scalar of the same type and kind as the output argument `x` representing
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and upper limit.<br>
    !>                                  Note that `gammaIncUpp` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>                                  Note that the procedure will abruptly end the program by calling `error stop`
    !>                                  if the computation of the Incomplete Gamma function fails to converge**.<br>
    !>
    !>  \interface{getGammaIncUppGil}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaGil, only: getGammaIncUppGil
    !>
    !>      gammaIncUpp = getGammaIncUppGil(x, kappa)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGammaIncLowGil](@ref pm_mathGammaGil::getGammaIncLowGil)<br>
    !>  [getGammaIncUppGil](@ref pm_mathGammaGil::getGammaIncUppGil)<br>
    !>  [setGammaIncGil](@ref pm_mathGammaGil::setGammaIncGil)<br>
    !>
    !>  \example{getGammaIncUppGil}
    !>  \include{lineno} example/pm_mathGammaGil/getGammaIncUppGil/main.F90
    !>  \compilef{getGammaIncUppGil}
    !>  \output{getGammaIncUppGil}
    !>  \include{lineno} example/pm_mathGammaGil/getGammaIncUppGil/main.out.F90
    !>  \postproc{getGammaIncUppGil}
    !>  \include{lineno} example/pm_mathGammaGil/getGammaIncUppGil/main.py
    !>  \vis{getGammaIncUppGil}
    !>  \image html pm_mathGammaGil/getGammaIncUppGil/getGammaIncUppGil.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaGil](@ref test_pm_mathGammaGil)
    !>
    !>  \final{getGammaIncUppGil}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncUppGil

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGammaIncUppGil_RK5(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppGil_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGammaIncUppGil_RK4(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppGil_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGammaIncUppGil_RK3(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppGil_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGammaIncUppGil_RK2(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppGil_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGammaIncUppGil_RK1(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppGil_RK1
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
    !>  Return the **regularized** Lower and Upper Incomplete Gamma function values
    !>  for the specified shape parameter (\f$\kappa\f$) and upper limit of the integral `x`.<br>
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
    !>  The Regularized Upper Incomplete Gamma Function can be readily obtained as,<br>
    !>
    !>  \f{equation}{
    !>      \large
    !>      Q(\kappa, x) = 1 - P(\kappa, x) ~,
    !>  \f}
    !>
    !>  However, the sake of numerical accuracy of either integral near `1`,
    !>  both lower and upper incomplete gamma values are returned by the procedures of this generic interface.<br>
    !>
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Lower Incomplete Gamma function also represents the Cumulative Distribution Function (CDF)
    !>  of the univariate Gamma distribution with the specified shape parameter and standardized `x` (with the scale parameter of unity).<br>
    !>
    !>  \param[out] gammaIncLow     :   The output scalar of the same type and kind as the input argument `x` representing
    !>                                  the Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncLow` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>  \param[out] gammaIncUpp     :   The output scalar of the same type and kind as the input argument `x` representing
    !>                                  the Upper Incomplete Gamma function for the specified `kappa` and lower limit.<br>
    !>                                  Note that `gammaIncUpp` is, by definition, always positive in the range \f$[0, 1]\f$.<br>
    !>  \param[in]  x               :   The input scalar of the type `real` of kind \RKALL,
    !>                                  representing the upper limit in the integral of the Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to zero if the algorithm succeeds to converge or a negative value if the algorithm fails to converge.<br>
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>
    !>  \interface{setGammaIncGil}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaGil, only: setGammaIncGil
    !>
    !>      call setGammaIncGil(gammaIncLow, gammaIncUpp, x, kappa, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < kappa` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGammaIncLowGil](@ref pm_mathGammaGil::getGammaIncLowGil)<br>
    !>  [setGammaIncGil](@ref pm_mathGammaGil::setGammaIncGil)<br>
    !>
    !>  \example{setGammaIncGil}
    !>  \include{lineno} example/pm_mathGammaGil/setGammaIncGil/main.F90
    !>  \compilef{setGammaIncGil}
    !>  \output{setGammaIncGil}
    !>  \include{lineno} example/pm_mathGammaGil/setGammaIncGil/main.out.F90
    !>  \postproc{setGammaIncGil}
    !>  \include{lineno} example/pm_mathGammaGil/setGammaIncGil/main.py
    !>  \vis
    !>  \image html pm_mathGammaGil/setGammaIncGil/setGammaIncGil.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaGil](@ref test_pm_mathGammaGil)
    !>
    !>
    !>  \final{setGammaIncGil}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncGil

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncGilDef_RK5(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncGilDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncGilDef_RK4(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncGilDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncGilDef_RK3(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncGilDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncGilDef_RK2(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncGilDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncLow, gammaIncUpp
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncGilDef_RK1(gammaIncLow, gammaIncUpp, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncGilDef_RK1
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

end module pm_mathGammaGil ! LCOV_EXCL_LINE
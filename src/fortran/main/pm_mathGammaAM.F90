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
!>          [setGammaIncLowAM](@ref pm_mathGammaAM::setGammaIncLowAM)
!>          to compute the Lower Incomplete Gamma function.<br>
!>      -#  **If ease of use matters** more than performance, use the function
!>          interfaces [getGammaIncLowAM](@ref pm_mathGammaAM::getGammaIncLowAM)
!>          to compute the Lower Incomplete Gamma function.<br>
!>
!>  \warning
!>  Although all generic interfaces of this module are available for  all processor `real` kinds,
!>  the accuracy and performance of the implemented algorithms are optimized for IEEE double precision.<br>
!>  In particular, the algorithms **may not accurately** compute the lower incomplete gamma function
!>  in extended precision (e.g., 128 bits) mode corresponding to \RKH kind type parameter.<br>
!>
!>  \note
!>  The computations of this module are explicitly based on the proposed approach by:<br>
!>  Abergel and Moisan, 2020, Algorithm 1006: Fast and Accurate Evaluation<br>
!>
!>  \see
!>  [pm_mathGamma](@ref pm_mathGamma) for detailed description of the (Regularized Incomplete) Gamma Function.<br>
!>
!>  \test
!>  [test_pm_mathGammaAM](@ref test_pm_mathGammaAM)
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

module pm_mathGammaAM

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathGammaAM"

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
    !>  \interface{getGammaIncLowAM}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaAM, only: getGammaIncLowAM
    !>
    !>      gammaIncLow = getGammaIncLowAM(x, kappa)
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
    !>  equivalents [setGammaIncLowAM](@ref setGammaIncLowAM).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [setGammaIncLowAM](@ref pm_mathGammaAM::setGammaIncLowAM)<br>
    !>  [getGammaIncLowAM](@ref pm_mathGammaAM::getGammaIncLowAM)<br>
    !>  [setGammaIncLowAM](@ref pm_mathGammaAM::setGammaIncLowAM)<br>
    !>
    !>  \example{getGammaIncLowAM}
    !>  \include{lineno} example/pm_mathGammaAM/getGammaIncLowAM/main.F90
    !>  \compilef{getGammaIncLowAM}
    !>  \output{getGammaIncLowAM}
    !>  \include{lineno} example/pm_mathGammaAM/getGammaIncLowAM/main.out.F90
    !>  \postproc{getGammaIncLowAM}
    !>  \include{lineno} example/pm_mathGammaAM/getGammaIncLowAM/main.py
    !>  \vis{getGammaIncLowAM}
    !>  \image html pm_mathGammaAM/getGammaIncLowAM/getGammaIncLowAM.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaAM](@ref test_pm_mathGammaAM)
    !>
    !>  \final{getGammaIncLowAM}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncLowAM

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncLowAM_RK5(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowAM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncLowAM_RK4(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowAM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncLowAM_RK3(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowAM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncLowAM_RK2(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowAM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncLow
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncLowAM_RK1(x, kappa) result(gammaIncLow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncLowAM_RK1
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
    !>  Note that this integral is bounded between zero and one (\f$[0,1]\f$).<br>
    !>  The regularized Upper Incomplete Gamma function also represents the *complement* of
    !>  the Cumulative Distribution Function (CDF) of the univariate Gamma distribution with
    !>  the specified shape parameter and standardized `x` (with the scale parameter of unity).
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
    !>  \interface{getGammaIncUppAM}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaAM, only: getGammaIncUppAM
    !>
    !>      gammaIncUpp = getGammaIncUppAM(x, kappa)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The procedures under this generic interface solely provide a more convenient method of calling the subroutine
    !>  equivalents [setGammaIncUppAM](@ref setGammaIncUppAM).<br>
    !>  As such they are slower than the corresponding subroutine version.<br>
    !>
    !>  \see
    !>  [setGammaIncUppAM](@ref pm_mathGammaAM::setGammaIncUppAM)<br>
    !>  [getGammaIncLowAM](@ref pm_mathGammaAM::getGammaIncLowAM)<br>
    !>  [setGammaIncLowAM](@ref pm_mathGammaAM::setGammaIncLowAM)<br>
    !>
    !>  \example{getGammaIncUppAM}
    !>  \include{lineno} example/pm_mathGammaAM/getGammaIncUppAM/main.F90
    !>  \compilef{getGammaIncUppAM}
    !>  \output{getGammaIncUppAM}
    !>  \include{lineno} example/pm_mathGammaAM/getGammaIncUppAM/main.out.F90
    !>  \postproc{getGammaIncUppAM}
    !>  \include{lineno} example/pm_mathGammaAM/getGammaIncUppAM/main.py
    !>  \vis{getGammaIncUppAM}
    !>  \image html pm_mathGammaAM/getGammaIncUppAM/getGammaIncUppAM.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaAM](@ref test_pm_mathGammaAM)
    !>
    !>  \remark
    !>  See Numerical Recipes by Press et al. 1992 for further details of the Incomplete Gamma function.
    !>
    !>  \final{getGammaIncUppAM}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface getGammaIncUppAM

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGammaIncUppAM_RK5(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppAM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGammaIncUppAM_RK4(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppAM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGammaIncUppAM_RK3(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppAM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGammaIncUppAM_RK2(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppAM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x, kappa
        real(RKG)                               :: gammaIncUpp
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGammaIncUppAM_RK1(x, kappa) result(gammaIncUpp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaIncUppAM_RK1
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
    !>  \param[in]  kappa           :   The input scalar of the same type and kind as `x`,
    !>                                  representing the shape parameter (\f$\kappa\f$) of the Lower Incomplete Gamma function \f$P(\kappa,x)\f$.
    !>  \param[out] info            :   The input scalar or array of the same shape as other input arguments of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to (positive) number of iterations taken for the algorithm to converge or its negative if the algorithm fails to converge.<br>
    !>                                  A convergence failure could happen if the input value for `kappa` is too large.<br>
    !>                                  **A negative value implies the lack of convergence.**<br>
    !>
    !>  \interface{setGammaIncLowAM}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaAM, only: setGammaIncLowAM
    !>
    !>      call setGammaIncLowAM(gammaIncLow, x, kappa, info)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the Gamma and related distributions.<br>
    !>  The logic behind pre-computing and passing `logGammaKappa = log_gamma(kappa)` is to speed up the calculations
    !>  since `log_gamma()` is computationally expensive and its recomputation can be avoided in repeated calls to
    !>  [setGammaIncLowAM](@ref pm_mathGammaAM::setGammaIncLowAM) with the same shape parameter but different `x`.
    !>
    !>  \see
    !>  [getGammaIncUppAM](@ref pm_mathGammaAM::getGammaIncUppAM)<br>
    !>  [setGammaIncUppAM](@ref pm_mathGammaAM::setGammaIncUppAM)<br>
    !>  [getGammaIncLowAM](@ref pm_mathGammaAM::getGammaIncLowAM)<br>
    !>
    !>  \example{setGammaIncLowAM}
    !>  \include{lineno} example/pm_mathGammaAM/setGammaIncLowAM/main.F90
    !>  \compilef{setGammaIncLowAM}
    !>  \output{setGammaIncLowAM}
    !>  \include{lineno} example/pm_mathGammaAM/setGammaIncLowAM/main.out.F90
    !>  \postproc{setGammaIncLowAM}
    !>  \include{lineno} example/pm_mathGammaAM/setGammaIncLowAM/main.py
    !>  \vis
    !>  \image html pm_mathGammaAM/setGammaIncLowAM/setGammaIncLowAM.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaAM](@ref test_pm_mathGammaAM)
    !>
    !>
    !>  \final{setGammaIncLowAM}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncLowAM

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncLowAM_RK5(gammaIncLow, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowAM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncLowAM_RK4(gammaIncLow, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowAM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncLowAM_RK3(gammaIncLow, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowAM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncLowAM_RK2(gammaIncLow, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowAM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncLowAM_RK1(gammaIncLow, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncLowAM_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncLow
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
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
    !>
    !>  \interface{setGammaIncUppAM}
    !>  \code{.F90}
    !>
    !>      use pm_mathGammaAM, only: setGammaIncUppAM
    !>
    !>      call setGammaIncUppAM(gammaIncUpp, x, logGammaKappa, kappa, info)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  These procedures are particularly useful and needed in computing the PDF of the Gamma and related distributions.<br>
    !>  The logic behind pre-computing and passing `logGammaKappa = log_gamma(kappa)` is to speed up the calculations
    !>  since `log_gamma()` is computationally expensive and its recomputation can be avoided in repeated calls to
    !>  [setGammaIncUppAM](@ref pm_mathGammaAM::setGammaIncUppAM) with the same shape parameter but different `x`.
    !>
    !>  \see
    !>  [getGammaIncUppAM](@ref pm_mathGammaAM::getGammaIncUppAM)<br>
    !>  [getGammaIncLowAM](@ref pm_mathGammaAM::getGammaIncLowAM)<br>
    !>  [setGammaIncLowAM](@ref pm_mathGammaAM::setGammaIncLowAM)<br>
    !>
    !>  \example{setGammaIncUppAM}
    !>  \include{lineno} example/pm_mathGammaAM/setGammaIncUppAM/main.F90
    !>  \compilef{setGammaIncUppAM}
    !>  \output{setGammaIncUppAM}
    !>  \include{lineno} example/pm_mathGammaAM/setGammaIncUppAM/main.out.F90
    !>  \postproc{setGammaIncUppAM}
    !>  \include{lineno} example/pm_mathGammaAM/setGammaIncUppAM/main.py
    !>  \vis{setGammaIncUppAM}
    !>  \image html pm_mathGammaAM/setGammaIncUppAM/setGammaIncUppAM.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathGammaAM](@ref test_pm_mathGammaAM)
    !>
    !>  \remark
    !>  See Numerical Recipes by Press et al. 1992 for further details of the Incomplete Gamma function.
    !>
    !>  \final{setGammaIncUppAM}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 12:36 pm, August 16, 2021, Dallas TX
    interface setGammaIncUppAM

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncUppAM_RK5(gammaIncUpp, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppAM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncUppAM_RK4(gammaIncUpp, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppAM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncUppAM_RK3(gammaIncUpp, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppAM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncUppAM_RK2(gammaIncUpp, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppAM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncUppAM_RK1(gammaIncUpp, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncUppAM_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaIncUpp
        real(RKG)   , intent(in)                :: x, logGammaKappa, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    ! This algorithm is taken from Abergel and Moisan, 2020, Algorithm 1006: Fast and Accurate Evaluation of a Generalized Incomplete Gamma Function
    ! It claims superb performance and accuracy. However, I could not substantiate their claims regarding accuracy, particularly,
    ! for very large `kappa` and `x` values `(kappa, x) > 1e7`.
    interface setGammaIncAM

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncAMDef_RK5(gammaInc, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncAMDef_RK4(gammaInc, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncAMDef_RK3(gammaInc, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncAMDef_RK2(gammaInc, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncAMDef_RK1(gammaInc, x, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGammaIncAMNXIK_RK5(gammaInc, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMNXIK_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, logGammaKappa
        integer(IK) , intent(in)                :: kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGammaIncAMNXIK_RK4(gammaInc, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMNXIK_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, logGammaKappa
        integer(IK) , intent(in)                :: kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGammaIncAMNXIK_RK3(gammaInc, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMNXIK_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, logGammaKappa
        integer(IK) , intent(in)                :: kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGammaIncAMNXIK_RK2(gammaInc, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMNXIK_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, logGammaKappa
        integer(IK) , intent(in)                :: kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGammaIncAMNXIK_RK1(gammaInc, x, logGammaKappa, kappa, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGammaIncAMNXIK_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)               :: gammaInc
        real(RKG)   , intent(in)                :: x, logGammaKappa
        integer(IK) , intent(in)                :: kappa
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !>  \endcond

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathGammaAM ! LCOV_EXCL_LINE
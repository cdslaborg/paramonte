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
!>  This module contains classes and procedures for computing various statistical quantities related to the mathematical <b>EggBox density function</b>.
!>
!>  \details
!>  The EggBox density function is frequently used in testing the efficiency of optimization and sampling algorithms.<br>
!>  Given a value \f$x\in\mathbb{R}^n\f$, the EggBox function with the <b>(intercept, location, scale, shape)</b> parameters \f$(\zeta, \mu, \sigma, \alpha)\f$ is defined as,
!>  \f{equation}{
!>      f(x) = \exp\left( \left[ \zeta + \prod_{i = 1}^{i = n} \cos\left(\pi\frac{x_i - \mu_i}{\sigma_i}\right)\right]^\alpha \right) ~,
!>  \f}
!>  where \f$\pi\f$ in the right-hand-side expression is the number \f$\ms{Pi}\f$.<br>
!>
!>  \test
!>  [test_pm_distEggBox](@ref test_pm_distEggBox)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distEggBox

    use pm_kind, only: IK, LK, RK, SK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distEggBox"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the EggBox density function at the specified input point `X`.
    !>
    !>  \details
    !>  See the documentation of [pm_distEggBox](@ref pm_distEggBox) for details of the EggBox density function.<br>
    !>
    !>  \param[in]  X               :   The input `contiguous` vector of size `(1:ndim)` of,
    !>                                  <ul>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ul>
    !>                                  containing the `ndim`-dimensional point at which the function must be evaluated.
    !>  \param[in]  mu              :   The input `contiguous` vector of the same type, kind, rank, and size as the input `X`,
    !>                                  representing the location parameter (\f$\mu\f$) of the density function.<br>
    !>                                  (**optional**, default = `0.`. It must be present **if and only if** the input argument `sigma` is also present.)
    !>  \param[in]  sigma           :   The input `contiguous` vector of the same type, kind, rank, and size as the input `X`,
    !>                                  representing the scale parameter (\f$\sigma\f$) of the density function.<br>
    !>                                  (**optional**, default = `1.` It must be present **if and only if** the input argument `mu` is also present.)
    !>  \param[in]  alpha           :   The input scalar of the same type and kind as the input `X`,
    !>                                  representing the shape parameter (\f$\alpha\f$) of the density function.<br>
    !>                                  (**optional**, default = `5.`)
    !>  \param[in]  zeta            :   The input scalar of the same type and kind as the input `X`,
    !>                                  representing the intercept parameter (\f$\zeta\f$) of the density function.<br>
    !>                                  (**optional**, default = `2.`)
    !>
    !>  \return
    !>  `logUDF`                    :   The output scalar of the same type and kind as the output argument `x`
    !>                                  representing the value of the density function at the specified location.
    !>
    !>  \interface{getEggBoxLogUDF}
    !>  \code{.F90}
    !>
    !>      use pm_distEggBox, only: getEggBoxLogUDF
    !>
    !>      logUDF = getEggBoxLogUDF(X(1:ndim), alpha = alpha, zeta = zeta)
    !>      logUDF = getEggBoxLogUDF(X(1:ndim), mu(1:ndim), sigma(1:ndim), alpha = alpha, zeta = zeta)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(X) == size(mu)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(X) == size(sigma)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \example{getEggBoxLogUDF}
    !>  \include{lineno} example/pm_distEggBox/getEggBoxLogUDF/main.F90
    !>  \compilef{getEggBoxLogUDF}
    !>  \output{getEggBoxLogUDF}
    !>  \include{lineno} example/pm_distEggBox/getEggBoxLogUDF/main.out.F90
    !>  \postproc{getEggBoxLogUDF}
    !>  \include{lineno} example/pm_distEggBox/getEggBoxLogUDF/main.py
    !>  \vis{getEggBoxLogUDF}
    !>  \image html pm_distEggBox/getEggBoxLogUDF/getEggBoxLogUDF.D1.RK.png width=700
    !>  \image html pm_distEggBox/getEggBoxLogUDF/getEggBoxLogUDF.D2.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distEggBox](@ref test_pm_distEggBox)
    !>
    !>  \finmain{getEggBoxLogUDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getEggBoxLogUDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getEggBoxLogUDFDDAZ_D1_RK5(X, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFDDAZ_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    , contiguous            :: X(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getEggBoxLogUDFDDAZ_D1_RK4(X, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFDDAZ_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    , contiguous            :: X(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getEggBoxLogUDFDDAZ_D1_RK3(X, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFDDAZ_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    , contiguous            :: X(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getEggBoxLogUDFDDAZ_D1_RK2(X, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFDDAZ_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    , contiguous            :: X(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getEggBoxLogUDFDDAZ_D1_RK1(X, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFDDAZ_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    , contiguous            :: X(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getEggBoxLogUDFMSAZ_D1_RK5(X, mu, sigma, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFMSAZ_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    , contiguous            :: X(:), mu(:), sigma(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getEggBoxLogUDFMSAZ_D1_RK4(X, mu, sigma, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFMSAZ_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    , contiguous            :: X(:), mu(:), sigma(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getEggBoxLogUDFMSAZ_D1_RK3(X, mu, sigma, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFMSAZ_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    , contiguous            :: X(:), mu(:), sigma(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getEggBoxLogUDFMSAZ_D1_RK2(X, mu, sigma, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFMSAZ_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    , contiguous            :: X(:), mu(:), sigma(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getEggBoxLogUDFMSAZ_D1_RK1(X, mu, sigma, alpha, zeta) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEggBoxLogUDFMSAZ_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    , contiguous            :: X(:), mu(:), sigma(:)
        real(RKC)   , intent(in)                , optional  :: alpha, zeta
        real(RKC)                                           :: logUDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distEggBox ! LCOV_EXCL_LINE
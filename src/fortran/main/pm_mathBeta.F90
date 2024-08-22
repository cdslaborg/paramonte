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
!>  This module contains classes and procedures for computing the mathematical Beta Function and its inverse.
!>
!>  \details
!>  The **Beta Function** is defined by the integral,
!>  \f{eqnarray}{
!>
!>      \mathrm{B}(\alpha, \beta)
!>      &=& \int_{0}^{1} ~ t^{\alpha-1} (1-t)^{\beta-1}~\mathrm{d}t \\
!>      &=& \frac{\Gamma(\alpha) \Gamma(\beta)}{\Gamma(\alpha + \beta)} ~,
!>
!>  \f}
!>  where \f$\alpha\f$ and \f$\beta\f$ are complex numbers such that \f$\Re \alpha > 0\f$ and  \f$\Re \beta > 0\f$.<br>
!>
!>  The **incomplete beta function** is a generalization of the beta function defined as,
!>  \f{eqnarray}{
!>
!>      \mathrm{B}(x; \alpha, \beta) = \int_{0}^{x} ~ t^{\alpha-1} (1-t)^{\beta-1}~\mathrm{d}t ~,~ x \in [0, 1]
!>
!>  \f}
!>  where for \f$x = 1\f$ becomes the same as the equation for the (complete) beta function,
!>  that is, \f$\mathrm{B}(x = 1;\alpha, \beta) = \mathrm{B}(\alpha, \beta)\f$.<br>
!>  Note that the relationship between the complete and the incomplete beta function is similar to
!>  the relationship between the complete and incomplete gamma function.<br>
!>
!>  The **regularized incomplete beta function** (or **Beta Function Ratio** or **regularized beta function** for short) is defined in terms of
!>  the incomplete beta function and the complete beta function as,
!>  \f{eqnarray}{
!>
!>      I_x(\alpha, \beta) = \frac{\mathrm{B}(x; \alpha, \beta)}{\mathrm{B}(\alpha, \beta)} ~.
!>
!>  \f}
!>  Note the regularized beta function has the limiting values,
!>  \f{eqnarray}{
!>
!>      I_0(\alpha, \beta) &=& 0 \\
!>      I_1(\alpha, \beta) &=& 1 ~.
!>
!>  \f}
!>
!>  The regularized incomplete beta function is the cumulative distribution function of the **beta distribution**,
!>  and is related to the cumulative distribution function of a random variable \f$X\f$ following a binomial
!>  distribution with probability of single success \f$p\f$ and number of Bernoulli trials \f$n\f$.
!>
!>  The **Regularized Inverse Incomplete Beta Function** \f$I_y(I_x; \alpha, \beta)\f$ is defined such that,
!>  \f{equation}{
!>      I_y\left(I_x(x; \alpha, \beta); \alpha, \beta\right) = x ~.
!>  \f}
!>  The **Regularized Inverse Incomplete Beta Function** is also the quantile function of the [Beta distribution](@ref pm_distBeta).<br>
!>
!>  \see
!>  [pm_mathGamma](@ref pm_mathGamma)<br>
!>  [pm_distBeta](@ref pm_distBeta)<br>
!>  [pm_distBeta](@ref pm_distBeta)<br>
!>  Newby, 1991, *The incomplete Beta Integral*, TH Eindhoven. THE/BDK/ORS, Vakgroep ORS: rapporten, Volume 9106.<br>
!>  Majumder and Bhattacharjee, 1973, *Algorithm AS 63: The incomplete Beta Integral*, Applied Statistics, Volume 22, Number 3, pages 409-411.<br>
!>
!>  \test
!>  [test_pm_mathBeta](@ref test_pm_mathBeta)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathBeta

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathBeta"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Beta Function \f$\mathrm{B}(\alpha, \beta)\f$
    !>  as defined in the details section of [pm_mathBeta](@ref pm_mathBeta).
    !>
    !>  \param[in]  alpha   :   The input scalar (or array of the same shape as other array-like arguments) of type `real` of kind \RKALL.<br>
    !>  \param[in]  beta    :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`.
    !>
    !>  \return
    !>  `logFuncBeta`       :   The output object of the same type, kind, and rank as highest-rank
    !>                          input argument containing the natural logarithm of the Beta Function.
    !>
    !>  \interface{getLogBeta}
    !>  \code{.F90}
    !>
    !>      use pm_mathBeta, only: getLogBeta
    !>
    !>      logFuncBeta = getLogBeta(alpha, beta)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 < alpha .and. 0 < beta` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [getBetaInc](@ref pm_mathBeta::getBetaInc)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{getLogBeta}
    !>  \include{lineno} example/pm_mathBeta/getLogBeta/main.F90
    !>  \compilef{getLogBeta}
    !>  \output{getLogBeta}
    !>  \include{lineno} example/pm_mathBeta/getLogBeta/main.out.F90
    !>  \postproc{getLogBeta}
    !>  \include{lineno} example/pm_mathBeta/getLogBeta/main.py
    !>  \vis{getLogBeta}
    !>  \image html pm_mathBeta/getLogBeta/getLogBeta.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathBeta](@ref test_pm_mathBeta)
    !>
    !>  \final{getLogBeta}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogBeta

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogBeta_RK5(alpha, beta) result(logFuncBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBeta_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: alpha, beta
        real(RKG)                   :: logFuncBeta
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogBeta_RK4(alpha, beta) result(logFuncBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBeta_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: alpha, beta
        real(RKG)                   :: logFuncBeta
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogBeta_RK3(alpha, beta) result(logFuncBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBeta_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: alpha, beta
        real(RKG)                   :: logFuncBeta
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogBeta_RK2(alpha, beta) result(logFuncBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBeta_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: alpha, beta
        real(RKG)                   :: logFuncBeta
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogBeta_RK1(alpha, beta) result(logFuncBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBeta_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: alpha, beta
        real(RKG)                   :: logFuncBeta
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **regularized** Incomplete Beta Function \f$I_x(\alpha, \beta)\f$ as defined in the details section of [pm_mathBeta](@ref pm_mathBeta).
    !>
    !>  \param[in]  x           :   The input scalar (or array of the same shape as other array-like arguments) of type `real` of kind \RKALL.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `x`.
    !>  \param[in]  beta        :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `x`.
    !>  \param[in]  signed      :   The input scalar (or array of the same shape as other array-like arguments) of type `logical` of default kind LK.<br>
    !>                              <ol>
    !>                                  <li>    If `signed = .false.`, the input `x` must be in range `0 <= x <= 1` and
    !>                                          the output `betaInc` will be the expected incomplete Beta function in range `0 <= betaInc <= 1`.<br>
    !>                                  <li>    If `signed = .true.`, then the following rules hold:
    !>                                          <ol>
    !>                                              <li>     If the condition `x < 0` holds, then the value `x = 1 - x < 0` will be used instead of `x`.<br>
    !>                                              <li>     If the output `betaInc` is near `1`, the output will be returned as `betaInc = betaInc - 1 < 0` instead of `betaInc`.<br>
    !>                                                      Therefore, the user is expected to be aware of the convention and apply the necessary correction (addition by `1`) before using the output value.<br>
    !>                                          </ol>
    !>                              </ol>
    !>                              Specifying `signed = .true.` can lead to considerably more accurate calculations (by orders of magnitudes) for values of `x` and `betaInc` that are near `1`.<br>
    !>                              The loss of precision near `1` occurs because of inadequate resolution of `real` number representations in digital computers near `1`
    !>                              which is orders of magnitude worse than the precision near `0`.<br>
    !>                              (**optional**, default = `.false.`, following the principle of least surprise.)
    !>  
    !>  \return 
    !>  `betaInc`               :   The output object of the same type, kind, and rank as highest-rank
    !>                              input argument containing the **regularized** Incomplete Beta Function.
    !>
    !>  \interface{getBetaInc}
    !>  \code{.F90}
    !>
    !>      use pm_mathBeta, only: getBetaInc
    !>
    !>      betaInc = getBetaInc(x, alpha, beta, signed = signed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 <= x .and. .not. signed .or. 0 <= x` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < alpha .and. 0 < beta` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The procedures under this generic interface will abort the program if the computation of
    !>  the continued fraction representation of the **regularized** beta function fails to converge.<br>
    !>  If error control is necessary, use the generic interface [seBetaInc](@ref pm_mathBeta::setBetaInc).<br>
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [getBetaInc](@ref pm_mathBeta::getBetaInc)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{getBetaInc}
    !>  \include{lineno} example/pm_mathBeta/getBetaInc/main.F90
    !>  \compilef{getBetaInc}
    !>  \output{getBetaInc}
    !>  \include{lineno} example/pm_mathBeta/getBetaInc/main.out.F90
    !>  \postproc{getBetaInc}
    !>  \include{lineno} example/pm_mathBeta/getBetaInc/main.py
    !>  \vis{getBetaInc}
    !>  \image html pm_mathBeta/getBetaInc/getBetaInc.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathBeta](@ref test_pm_mathBeta)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBetaInc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getBetaInc_RK5(x, alpha, beta, signed) result(betaInc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInc_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                               :: betaInc
        real(RKG)   , intent(in)                :: x, alpha, beta
        logical(LK) , intent(in)    , optional  :: signed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getBetaInc_RK4(x, alpha, beta, signed) result(betaInc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInc_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                               :: betaInc
        real(RKG)   , intent(in)                :: x, alpha, beta
        logical(LK) , intent(in)    , optional  :: signed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getBetaInc_RK3(x, alpha, beta, signed) result(betaInc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInc_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                               :: betaInc
        real(RKG)   , intent(in)                :: x, alpha, beta
        logical(LK) , intent(in)    , optional  :: signed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getBetaInc_RK2(x, alpha, beta, signed) result(betaInc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInc_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                               :: betaInc
        real(RKG)   , intent(in)                :: x, alpha, beta
        logical(LK) , intent(in)    , optional  :: signed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getBetaInc_RK1(x, alpha, beta, signed) result(betaInc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInc_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                               :: betaInc
        real(RKG)   , intent(in)                :: x, alpha, beta
        logical(LK) , intent(in)    , optional  :: signed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **regularized** Incomplete Beta Function \f$I_x(\alpha, \beta)\f$ as defined in the details section of [pm_mathBeta](@ref pm_mathBeta).
    !>
    !>  \param[in]  betaInc     :   The output scalar (or array of the same shape as other array-like arguments)
    !>                              of the same type and kind as the input argument `x` containing the **regularized** Incomplete Beta Function.
    !>  \param[in]  x           :   The input scalar (or array of the same shape as other array-like arguments) of,<br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the value at which the function must be computed.
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `x`.
    !>  \param[in]  beta        :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `x`.
    !>  \param[in]  logFuncBeta :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `x`,
    !>                              representing the natural logarithm of the Beta Function as returned by [getLogBeta(alpha, beta)](@ref pm_mathBeta::getLogBeta).<br>
    !>                              Providing this argument can lead to significant runtime performance gains when the interface
    !>                              is called repeatedly for different `x` values but with identical `alpha` and `beta`.
    !>  \param[in]  signed      :   The input scalar (or array of the same shape as other array-like arguments) of type `logical` of default kind LK.<br>
    !>                              <ol>
    !>                                  <li>    If `signed = .false.`, the input `x` must be in range `0 <= x <= 1` and
    !>                                          the output `betaInc` will be the expected incomplete Beta function in range `0 <= betaInc <= 1`.<br>
    !>                                  <li>    If `signed = .true.`, then the following rules hold:
    !>                                          <ol>
    !>                                              <li>     If the condition `x < 0` holds, then the value `x = 1 - x < 0` will be used instead of `x`.<br>
    !>                                              <li>     If the output `betaInc` is near `1`, the output will be returned as `betaInc = betaInc - 1 < 0` instead of `betaInc`.<br>
    !>                                                      Therefore, the user is expected to be aware of the convention and apply the necessary correction (addition by `1`) before using the output value.<br>
    !>                                          </ol>
    !>                              </ol>
    !>                              Specifying `signed = .true.` can lead to considerably more accurate calculations (by orders of magnitudes) for values of `x` and `betaInc` that are near `1`.<br>
    !>                              The loss of precision near `1` occurs because of inadequate resolution of `real` number representations in digital computers near `1`
    !>                              which is orders of magnitude worse than the precision near `0`.<br>
    !>  \param[in]  reltol      :   The input scalar argument of the same type and kind as `betaInc`, representing the relative tolerance of integration in the definition of the Incomplete Beta function.<br>
    !>                              If the relative accuracy of integration reaches a value below this threshold, the computation of the Incomplete Beta function is assumed to have converged.<br>
    !>                              Specifying this argument will lead to brute-force numerical integration of the integrand of the Beta function to compute the integral.<br>
    !>                              This approach can be useful when the continued-fraction representation of the Incomplete Beta Function fails to converge or is extremely costly or a bound on the numerical error is needed.<br>
    !>                              Convergence failure or costly computations frequently occur for large values of `alpha` and `beta`.<br>
    !>                              This argument can be set to any positive value significantly larger (e.g., \f$\times10000\f$) than `epsilon(real(0., kind(reltol))`.<br>
    !>                              A good rule of thumb is to set `reltol = epsilon(real(0, kind(integral)))**(2./3.)`.<br>
    !>                              The integration result is generally orders of magnitude more precise than the specified `reltol`.<br>
    !>                              (**optional**. If missing, the continued-fraction representation of the Incomplete Beta function will be used to approximate the output.)
    !>  \param[out] info        :   The output scalar (or array of the same shape as other array-like arguments) of type `integer` of default kind \IK.<br>
    !>                              It is `0` <b>if and only if</b> the computation converges, otherwise,<br>
    !>                              <ol>
    !>                                  <li>    if the input argument `reltol` is missing, it is set to `1` to signal lack of convergence in the computation of the continued-fraction representation of the Incomplete Beta function.<br>
    !>                                          Convergence fails if `alpha` or `beta` are too large, in which case a brute-force integration method by specifying the input argument `reltol` might resolve it.<br>
    !>                                  <li>    if the input argument `reltol` is present, it is set to the non-zero `info` output of [getQuadErr](@ref pm_quadPack::getQuadErr) to signal lack of convergence in the computation of the integral.<br>
    !>                                          See the documentation of the output argument `info` of [getQuadErr](@ref pm_quadPack::getQuadErr) for possible output error codes and their meanings.<br>
    !>                              </ol>
    !>                              The failure may be resolved by increasing the preset value of `MAX_ITER` in the implementation.<br>
    !>
    !>  \return
    !>
    !>  \interface{setBetaInc}
    !>  \code{.F90}
    !>
    !>      use pm_mathBeta, only: setBetaInc
    !>
    !>      call setBetaInc(betaInc, x, alpha, beta, logFuncBeta, signed, info)
    !>      call setBetaInc(betaInc, x, alpha, beta, logFuncBeta, reltol, signed, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 <= x .and. x <= 1` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < alpha .and. 0 < beta` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the input argument `reltol` is present.<br>
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [setBetaInc](@ref pm_mathBeta::setBetaInc)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{setBetaInc}
    !>  \include{lineno} example/pm_mathBeta/setBetaInc/main.F90
    !>  \compilef{setBetaInc}
    !>  \output{setBetaInc}
    !>  \include{lineno} example/pm_mathBeta/setBetaInc/main.out.F90
    !>  \postproc{setBetaInc}
    !>  \include{lineno} example/pm_mathBeta/setBetaInc/main.py
    !>  \vis{setBetaInc}
    !>  \image html pm_mathBeta/setBetaInc/setBetaInc.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathBeta](@ref test_pm_mathBeta)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBetaInc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBetaIncDef_RK5(betaInc, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBetaIncDef_RK4(betaInc, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBetaIncDef_RK3(betaInc, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBetaIncDef_RK2(betaInc, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBetaIncDef_RK1(betaInc, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBetaIncGK21_RK5(betaInc, x, alpha, beta, logFuncBeta, reltol, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncGK21_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        real(RKG)       , intent(in)                :: reltol
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBetaIncGK21_RK4(betaInc, x, alpha, beta, logFuncBeta, reltol, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncGK21_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        real(RKG)       , intent(in)                :: reltol
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBetaIncGK21_RK3(betaInc, x, alpha, beta, logFuncBeta, reltol, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncGK21_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        real(RKG)       , intent(in)                :: reltol
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBetaIncGK21_RK2(betaInc, x, alpha, beta, logFuncBeta, reltol, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncGK21_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        real(RKG)       , intent(in)                :: reltol
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBetaIncGK21_RK1(betaInc, x, alpha, beta, logFuncBeta, reltol, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaIncGK21_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)       , intent(out)               :: betaInc
        real(RKG)       , intent(in)                :: logFuncBeta
        real(RKG)       , value                     :: x, alpha, beta
        real(RKG)       , intent(in)                :: reltol
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **regularized** Inverse Incomplete Beta Function \f$I_x(\alpha, \beta)\f$ as defined in the details section of [pm_mathBeta](@ref pm_mathBeta).
    !>
    !>  \param[in]  betaInc     :   The input scalar (or array of the same shape as other array-like arguments) of type `real` of kind \RKALL.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `betaInc`.
    !>  \param[in]  beta        :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `betaInc`.
    !>  \param[in]  signed      :   The input scalar (or array of the same shape as other array-like arguments) of type `logical` of default kind LK.<br>
    !>                              <ol>
    !>                                  <li>    If `signed = .false.`, the input `betaInc` must be in range `0 <= betaInc <= 1` and
    !>                                          the output `betaInv` will be the expected incomplete Beta function in range `0 <= betaInv <= 1`.<br>
    !>                                  <li>    If `signed = .true.`, then the following rules hold:
    !>                                          <ol>
    !>                                              <li>     If the condition `betaInc < 0` holds, then the value `betaInc = 1 - betaInc < 0` will be used instead of `betaInc`.<br>
    !>                                              <li>     If the output `betaInv` is near `1`, the output will be returned as `betaInv = betaInv - 1 < 0` instead of `betaInv`.<br>
    !>                                                      Therefore, the user is expected to be aware of the convention and apply the necessary correction (addition by `1`) before using the output value.<br>
    !>                                          </ol>
    !>                              </ol>
    !>                              Specifying `signed = .true.` can lead to considerably more accurate calculations (by orders of magnitudes) for values of `betaInc` and `betaInv` that are near `1`.<br>
    !>                              The loss of precision near `1` occurs because of inadequate resolution of `real` number representations in digital computers near `1`
    !>                              which is orders of magnitude worse than the precision near `0`.<br>
    !>                              (**optional**, default = `.false.`, following the principle of least surprise.)
    !>
    !>  \return
    !>  `betaInv`               :   The output object of the same type, kind, and rank as highest-rank
    !>                              input argument containing the **regularized** Inverse Incomplete Beta Function.
    !>
    !>  \interface{getBetaInv}
    !>  \code{.F90}
    !>
    !>      use pm_mathBeta, only: getBetaInv
    !>
    !>      betaInv = getBetaInv(betaInc, alpha, beta, signed = signed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 <= betaInc .and. 0 <= betaInc` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < alpha .and. 0 < beta` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The procedures under this generic interface will abort the program if the computation of
    !>  the continued fraction representation of the **regularized** beta function fails to converge.
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [getBetaInv](@ref pm_mathBeta::getBetaInv)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{getBetaInv}
    !>  \include{lineno} example/pm_mathBeta/getBetaInv/main.F90
    !>  \compilef{getBetaInv}
    !>  \output{getBetaInv}
    !>  \include{lineno} example/pm_mathBeta/getBetaInv/main.out.F90
    !>  \postproc{getBetaInv}
    !>  \include{lineno} example/pm_mathBeta/getBetaInv/main.py
    !>  \vis{getBetaInv}
    !>  \image html pm_mathBeta/getBetaInv/getBetaInv.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathBeta](@ref test_pm_mathBeta)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBetaInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBetaInv_RK5(betaInc, alpha, beta, signed) result(betaInv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInv_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                                   :: betaInv
        real(RKG)       , intent(in)                :: betaInc, alpha, beta
        logical(LK)     , intent(in)    , optional  :: signed
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBetaInv_RK4(betaInc, alpha, beta, signed) result(betaInv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInv_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                                   :: betaInv
        real(RKG)       , intent(in)                :: betaInc, alpha, beta
        logical(LK)     , intent(in)    , optional  :: signed
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBetaInv_RK3(betaInc, alpha, beta, signed) result(betaInv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInv_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                                   :: betaInv
        real(RKG)       , intent(in)                :: betaInc, alpha, beta
        logical(LK)     , intent(in)    , optional  :: signed
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBetaInv_RK2(betaInc, alpha, beta, signed) result(betaInv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInv_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                                   :: betaInv
        real(RKG)       , intent(in)                :: betaInc, alpha, beta
        logical(LK)     , intent(in)    , optional  :: signed
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBetaInv_RK1(betaInc, alpha, beta, signed) result(betaInv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaInv_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                                   :: betaInv
        real(RKG)       , intent(in)                :: betaInc, alpha, beta
        logical(LK)     , intent(in)    , optional  :: signed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **regularized** Inverse Incomplete Beta Function \f$I_x(\alpha, \beta)\f$ as defined in the details section of [pm_mathBeta](@ref pm_mathBeta).
    !>
    !>  \param[in]  betaInv     :   The output scalar (or array of the same shape as other array-like arguments)
    !>                              of the same type and kind as the input argument `betaInc` containing the **regularized** Inverse Incomplete Beta Function.
    !>  \param[in]  betaInc     :   The input scalar (or array of the same shape as other array-like arguments) of,<br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the value at which the function must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `betaInc`.<br>
    !>  \param[in]  beta        :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `betaInc`.<br>
    !>  \param[in]  logFuncBeta :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `betaInc`,
    !>                              representing the natural logarithm of the Beta Function as returned by [getLogBeta(alpha, beta)](@ref pm_mathBeta::getLogBeta).<br>
    !>                              Providing this argument can lead to significant runtime performance gains when the interface
    !>                              is called repeatedly for different `betaInc` values but with identical `alpha` and `beta`.
    !>  \param[in]  signed      :   The input scalar (or array of the same shape as other array-like arguments) of type `logical` of default kind LK.<br>
    !>                              <ol>
    !>                                  <li>    If `signed = .false.`, the input `betaInc` must be in range `0 <= betaInc <= 1` and
    !>                                          the output `betaInv` will be the expected incomplete Beta function in range `0 <= betaInv <= 1`.<br>
    !>                                  <li>    If `signed = .true.`, then the following rules hold:
    !>                                          <ol>
    !>                                              <li>     If the condition `betaInc < 0` holds, then the value `betaInc = 1 - betaInc < 0` will be used instead of `betaInc`.<br>
    !>                                              <li>     If the output `betaInv` is near `1`, the output will be returned as `betaInv = betaInv - 1 < 0` instead of `betaInv`.<br>
    !>                                                      Therefore, the user is expected to be aware of the convention and apply the necessary correction (addition by `1`) before using the output value.<br>
    !>                                          </ol>
    !>                              </ol>
    !>                              Specifying `signed = .true.` can lead to considerably more accurate calculations (by orders of magnitudes) for values of `betaInc` and `betaInv` that are near `1`.<br>
    !>                              The loss of precision near `1` occurs because of inadequate resolution of `real` number representations in digital computers near `1`
    !>                              which is orders of magnitude worse than the precision near `0`.<br>
    !>  \param[out] info        :   The output scalar (or array of the same shape as other array-like arguments) of type `integer` of default kind \IK.<br>
    !>                              <ol>
    !>                                  <li>    It is set to `0` <b>if and only if</b> the algorithm returns successfully.<br>
    !>                                  <li>    It is set to a non-zero <b>if and only if</b> the computation of the inverse Incomplete Beta Function at the specified value fails to converge.<br>
    !>                                          See the documentation of the corresponding argument `info` of [setBetaInc](@ref pm_mathBeta::setBetaInc) for possible output values and their meanings.<br>
    !>                              </ol>
    !>
    !>  \return
    !>
    !>  \interface{setBetaInv}
    !>  \code{.F90}
    !>
    !>      use pm_mathBeta, only: setBetaInv
    !>
    !>      call setBetaInv(betaInv, betaInc, alpha, beta, logFuncBeta, signed, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 <= betaInc .and. betaInc <= 1` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < alpha .and. 0 < beta` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [setBetaInv](@ref pm_mathBeta::setBetaInv)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{setBetaInv}
    !>  \include{lineno} example/pm_mathBeta/setBetaInv/main.F90
    !>  \compilef{setBetaInv}
    !>  \output{setBetaInv}
    !>  \include{lineno} example/pm_mathBeta/setBetaInv/main.out.F90
    !>  \postproc{setBetaInv}
    !>  \include{lineno} example/pm_mathBeta/setBetaInv/main.py
    !>  \vis{setBetaInv}
    !>  \image html pm_mathBeta/setBetaInv/setBetaInv.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathBeta](@ref test_pm_mathBeta)
    !>
    !>  \final{setBetaInv}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBetaInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBetaInv_RK5(betaInv, betaInc, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaInv_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)       , intent(out)               :: betaInv
        real(RKG)       , value                     :: betaInc, alpha, beta, logFuncBeta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBetaInv_RK4(betaInv, betaInc, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaInv_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)       , intent(out)               :: betaInv
        real(RKG)       , value                     :: betaInc, alpha, beta, logFuncBeta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBetaInv_RK3(betaInv, betaInc, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaInv_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)       , intent(out)               :: betaInv
        real(RKG)       , value                     :: betaInc, alpha, beta, logFuncBeta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBetaInv_RK2(betaInv, betaInc, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaInv_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)       , intent(out)               :: betaInv
        real(RKG)       , value                     :: betaInc, alpha, beta, logFuncBeta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBetaInv_RK1(betaInv, betaInc, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaInv_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)       , intent(out)               :: betaInv
        real(RKG)       , value                     :: betaInc, alpha, beta, logFuncBeta
        logical(LK)     , intent(in)                :: signed
        integer(IK)     , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathBeta
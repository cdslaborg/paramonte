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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>(Truncated) Power distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>(Truncated) Power distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function **(ICDF)** or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the <b>(Truncated) Power distribution</b> over a strictly-positive support \f$x \in [x_\mathrm{min}, x_\mathrm{max}]\f$
!>  is defined with the three <b>(shape, scale, scale)</b> parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ as,<br>
!>  \f{equation}{
!>      \large
!>      \pi(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~ x^{\alpha - 1} ~,
!>  \f}
!>  where \f$\mathbf{0 < \alpha < +\infty}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\eta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPower::getPowerLogPDFNF) of the PDF,<br>
!>  \f{equation}{
!>      \large
!>      \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~,
!>  \f}
!>
!>  When \f$x_\mathrm{min} \rightarrow 0\f$, the <b>Truncated Power distribution</b> simplifies to the **Power Distribution** with PDF.<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{min} \rightarrow 0} \pi(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{\alpha}{x_\mathrm{max}^\alpha} x^{\alpha - 1} ~,
!>  \f}
!>
!>  The equation for \f$\eta(\cdot)\f$ for the Power distribution simplifies to,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{min} \rightarrow 0} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \alpha x_\mathrm{max}^{-\alpha} ~,~ 0 < \alpha < +\infty.
!>  \f}
!>
!>  The corresponding **CDF** of the (Truncated) Power distribution is given by,<br>
!>  \f{equation}{
!>      \large
!>      \mathrm{CDF}(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) \bigg[\bigg(\frac{x}{x_\mathrm{min}}\bigg)^\alpha - 1\bigg] ~,
!>  \f}
!>  where \f$\mathbf{0 < \alpha < +\infty}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPower::getPowerLogCDFNF) of the CDF,<br>
!>  \f{eqnarray}{
!>      \large
!>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})
!>      &=& \frac{x_\mathrm{min}^\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~, \\
!>      &=& \frac{x_\mathrm{min}^\alpha}{\alpha} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~,
!>  \f}
!>  where \f$\eta(\cdot)\f$ is the [normalization factor](@ref pm_distPower::getPowerLogPDFNF) of the PDF.<br>
!>
!>  When \f$x_\mathrm{min} \rightarrow 0\f$, the <b>Truncated Power distribution</b> simplifies to the **Power Distribution** with CDF,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{min} \rightarrow 0} \mathrm{CDF}(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \bigg(\frac{x}{x_\mathrm{max}}\bigg)^\alpha ~,
!>  \f}
!>  with,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{min} \rightarrow 0} \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = x_\mathrm{max}^{-\alpha} ~.
!>  \f}
!>
!>  The corresponding **Inverse CDF** of the (Truncated) Power distribution is given by,
!>  \f{equation}{
!>      \large
!>      Q(\mathrm{CDF}(x); \alpha, x_\mathrm{min}, x_\mathrm{max}) \equiv x = x_\mathrm{min} \bigg(1 + \frac{\mathrm{CDF}(x)}{\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})}\bigg)^{\frac{1}{\alpha}} ~,
!>  \f}
!>  where \f$\mathbf{0 < \alpha < +\infty}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq Q(\mathrm{CDF}(x); \alpha, x_\mathrm{min}, x_\mathrm{max}) \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPower::getPowerLogCDFNF) of the CDF,<br>
!>  \f{equation}{
!>      \large
!>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{x_\mathrm{min}^\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~.
!>  \f}
!>
!>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the **Truncated Power distribution** simplifies to the **Power Distribution**, with,<br>
!>  \f{equation}{
!>      \large
!>      Q(\mathrm{CDF}(x); \alpha, x_\mathrm{min}, x_\mathrm{max}) \equiv x = x_\mathrm{max} \big(\mathrm{CDF}(x)\big)^{\frac{1}{\alpha}} = \bigg(\frac{\mathrm{CDF}(x)}{\zeta(\alpha, x_\mathrm{max})}\bigg)^{\frac{1}{\alpha}} ~,
!>  \f}
!>
!>  **Random Number Generation**<br>
!>
!>  Assuming that \f$U \in [0, 1)\f$ is a uniformly-distributed random variate, the transformed random variable,<br>
!>  \f{equation}{
!>      \large
!>      x = x_\mathrm{min} \bigg(1 + \frac{U}{\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})}\bigg)^{\frac{1}{\alpha}} ~,
!>  \f}
!>  follows a **Truncated Power distribution** with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max}\f$)
!>  where \f$\mathbf{0 < \alpha < +\infty}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPower::getPowerLogCDFNF) of the CDF,<br>
!>  \f{equation}{
!>      \large
!>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{x_\mathrm{min}^\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~.
!>  \f}
!>
!>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the **Truncated Power distribution** simplifies to the **Power Distribution**, with,<br>
!>  \f{equation}{
!>      \large
!>      x = x_\mathrm{max} U^{\frac{1}{\alpha}} = \bigg(\frac{U}{\zeta(\alpha, x_\mathrm{max})}\bigg)^{\frac{1}{\alpha}} ~,
!>  \f}
!>
!>  \remark
!>  See [pm_distPareto](@ref pm_distPareto) for the case of \f$-\infty < \alpha < 0\f$.<br>
!>  See [pm_distPoweto](@ref pm_distPoweto) for the case of \f$\alpha = 0\f$.<br>
!>
!>  \see
!>  [pm_distPower](@ref pm_distPower)<br>
!>  [pm_distPareto](@ref pm_distPareto)<br>
!>  [pm_distPoweto](@ref pm_distPoweto)<br>
!>
!>  \test
!>  [test_pm_distPower](@ref test_pm_distPower)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distPower

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distPower"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the PDF
    !>  of the (Truncated) Power distribution for an input parameter set \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  The normalization factor \f$\eta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ of the PDF of the Power distribution is,
    !>
    !>  \f{equation}{
    !>
    !>      \large
    !>      \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~,~ 0 < \alpha < +\infty
    !>
    !>  \f}
    !>
    !>  When \f$x_\mathrm{min} \rightarrow 0\f$, the <b>Truncated Power distribution</b> simplifies to the **Power Distribution**.<br>
    !>  The equation for \f$\eta(\cdot)\f$ for the Power distribution simplifies to,<br>
    !>
    !>  \f{equation}{
    !>      \large
    !>      \lim_{x_\mathrm{min} \rightarrow 0} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \alpha x_\mathrm{max}^{-\alpha} ~,~ 0 < \alpha < +\infty.
    !>  \f}
    !>
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factor of the PDF of the (Truncated) Power distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the properties of the (Truncated) Power distribution to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha   :   The input scalar or array of the same shape as other array-like arguments, of type `real` of kind \RKALL,
    !>                          containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                          (**optional**, default = \f$-\infty\f$)
    !>  \param[in]  logMaxX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logPDFNF`          :   The output scalar or array of the same shape as any input array-like argument,
    !>                          of the same type and kind as the input argument `alpha`, containing the natural logarithm
    !>                          of the normalization factor of the PDF of (Truncated) Power distribution.<br>
    !>
    !>  \interface{getPowerLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: getPowerLogPDFNF
    !>
    !>      logPDFNF = getPowerLogPDFNF(alpha, logMaxX) ! Power distribution.
    !>      logPDFNF = getPowerLogPDFNF(alpha, logMinX, logMaxX) ! Truncated Power distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha` must hold for the corresponding input arguments.<br>
    !>  The condition `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowerLogPDF](@ref pm_distPower::getPowerLogPDF)<br>
    !>  [setPowerLogPDF](@ref pm_distPower::setPowerLogPDF)<br>
    !>
    !>  \example{getPowerLogPDFNF}
    !>  \include{lineno} example/pm_distPower/getPowerLogPDFNF/main.F90
    !>  \compilef{getPowerLogPDFNF}
    !>  \output{getPowerLogPDFNF}
    !>  \include{lineno} example/pm_distPower/getPowerLogPDFNF/main.out.F90
    !>  \postproc{getPowerLogPDFNF}
    !>  \include{lineno} example/pm_distPower/getPowerLogPDFNF/main.py
    !>  \vis{getPowerLogPDFNF}
    !>  \image html pm_distPower/getPowerLogPDFNF/getPowerLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowerLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowerLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogPDFNFALD_RK5(alpha, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogPDFNFALD_RK4(alpha, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogPDFNFALD_RK3(alpha, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogPDFNFALD_RK2(alpha, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogPDFNFALD_RK1(alpha, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogPDFNFALL_RK5(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogPDFNFALL_RK4(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogPDFNFALL_RK3(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogPDFNFALL_RK2(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogPDFNFALL_RK1(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFNFALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the
    !>  (Truncated) Power distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.
    !>
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$)
    !>  \param[in]  logMaxX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `logx`, containing the natural logarithm of the PDF of the distribution.<br>
    !>
    !>  \interface{getPowerLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: getPowerLogPDF
    !>
    !>      logPDF = getPowerLogPDF(logx, alpha, logMaxX) ! Power distribution.
    !>      logPDF = getPowerLogPDF(logx, alpha, logMinX, logMaxX) ! Truncated Power distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < alpha` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  The conditions `logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowerLogPDF](@ref pm_distPower::setPowerLogPDF)<br>
    !>
    !>  \example{getPowerLogPDF}
    !>  \include{lineno} example/pm_distPower/getPowerLogPDF/main.F90
    !>  \compilef{getPowerLogPDF}
    !>  \output{getPowerLogPDF}
    !>  \include{lineno} example/pm_distPower/getPowerLogPDF/main.out.F90
    !>  \postproc{getPowerLogPDF}
    !>  \include{lineno} example/pm_distPower/getPowerLogPDF/main.py
    !>  \vis{getPowerLogPDF}
    !>  \image html pm_distPower/getPowerLogPDF/getPowerLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowerLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowerLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogPDFALD_RK5(logx, alpha, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogPDFALD_RK4(logx, alpha, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogPDFALD_RK3(logx, alpha, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogPDFALD_RK2(logx, alpha, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogPDFALD_RK1(logx, alpha, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogPDFALL_RK5(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogPDFALL_RK4(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogPDFALL_RK3(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogPDFALL_RK2(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogPDFALL_RK1(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogPDFALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the
    !>  (Truncated) Power distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logPDFNF    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor of the PDF of the (Truncated) Power distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getPowerLogPDFNF(alpha, logMinX, logMaxX)](@ref pm_distPower::getPowerLogPDFNF).<br>
    !>
    !>  \interface{setPowerLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: setPowerLogPDF
    !>
    !>      call setPowerLogPDF(logPDF, logx, alpha, logPDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Be warned that the procedures under this generic interface have no mechanism for checking for the consistency of the specified input arguments with each other.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowerLogPDF](@ref pm_distPower::setPowerLogPDF)<br>
    !>
    !>  \example{setPowerLogPDF}
    !>  \include{lineno} example/pm_distPower/setPowerLogPDF/main.F90
    !>  \compilef{setPowerLogPDF}
    !>  \output{setPowerLogPDF}
    !>  \include{lineno} example/pm_distPower/setPowerLogPDF/main.out.F90
    !>  \postproc{setPowerLogPDF}
    !>  \include{lineno} example/pm_distPower/setPowerLogPDF/main.py
    !>  \vis{setPowerLogPDF}
    !>  \image html pm_distPower/setPowerLogPDF/setPowerLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setPowerLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowerLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogPDF_RK5(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogPDF_RK4(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogPDF_RK3(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogPDF_RK2(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogPDF_RK1(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the CDF
    !>  of the (Truncated) Power distribution for an input parameter set \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  The normalization factor \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ of the CDF of the Truncated Power distribution is,
    !>
    !>  \f{eqnarray}{
    !>      \large
    !>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})
    !>      &=& \frac{x_\mathrm{min}^\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~, \\
    !>      &=& \frac{x_\mathrm{min}^\alpha}{\alpha} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~,
    !>  \f}
    !>
    !>  where \f$\eta(\cdot)\f$ is the [normalization factor](@ref pm_distPower::getPowerLogPDFNF) of the PDF.<br>
    !>  The equation for \f$\zeta(\cdot)\f$ for the Power distribution simplifies to,<br>
    !>
    !>  \f{equation}{
    !>      \large
    !>      \lim_{x_\mathrm{min} \rightarrow 0} \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = x_\mathrm{max}^{-\alpha} ~.
    !>  \f}
    !>
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factor of the CDF of the (Truncated) Power distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the properties of the (Truncated) Power distribution to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha   :   The input scalar or array of the same shape as other array-like arguments, of type `real` of kind \RKALL,
    !>                          containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                          (**optional**, default = \f$-\infty\f$)
    !>  \param[in]  logMaxX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logCDFNF`          :   The output scalar or array of the same shape as any input array-like argument,
    !>                          of the same type and kind as the input argument `alpha`, containing the natural logarithm
    !>                          of the normalization factor of the CDF of (Truncated) Power distribution.<br>
    !>
    !>  \interface{getPowerLogCDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: getPowerLogCDFNF
    !>
    !>      logCDFNF = getPowerLogCDFNF(alpha, logMaxX) ! Power distribution.
    !>      logCDFNF = getPowerLogCDFNF(alpha, logMinX, logMaxX) ! Truncated Power distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha` must hold for the corresponding input arguments.<br>
    !>  The condition `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowerLogCDF](@ref pm_distPower::getPowerLogCDF)<br>
    !>  [setPowerLogCDF](@ref pm_distPower::setPowerLogCDF)<br>
    !>
    !>  \example{getPowerLogCDFNF}
    !>  \include{lineno} example/pm_distPower/getPowerLogCDFNF/main.F90
    !>  \compilef{getPowerLogCDFNF}
    !>  \output{getPowerLogCDFNF}
    !>  \include{lineno} example/pm_distPower/getPowerLogCDFNF/main.out.F90
    !>  \postproc{getPowerLogCDFNF}
    !>  \include{lineno} example/pm_distPower/getPowerLogCDFNF/main.py
    !>  \vis{getPowerLogCDFNF}
    !>  \image html pm_distPower/getPowerLogCDFNF/getPowerLogCDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowerLogCDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowerLogCDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogCDFNFALD_RK5(alpha, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogCDFNFALD_RK4(alpha, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogCDFNFALD_RK3(alpha, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogCDFNFALD_RK2(alpha, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogCDFNFALD_RK1(alpha, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogCDFNFALL_RK5(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogCDFNFALL_RK4(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogCDFNFALL_RK3(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogCDFNFALL_RK2(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogCDFNFALL_RK1(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFNFALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Cumulative Distribution Function (CDF) of the
    !>  (Truncated) Power distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.
    !>
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$)
    !>  \param[in]  logMaxX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logCDF`                :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `logx`, containing the natural logarithm of the CDF of the distribution.<br>
    !>
    !>  \interface{getPowerLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: getPowerLogCDF
    !>
    !>      logCDF = getPowerLogCDF(logx, alpha, logMaxX) ! Power distribution.
    !>      logCDF = getPowerLogCDF(logx, alpha, logMinX, logMaxX) ! Truncated Power distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < alpha` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  The conditions `logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowerLogCDF](@ref pm_distPower::setPowerLogCDF)<br>
    !>
    !>  \example{getPowerLogCDF}
    !>  \include{lineno} example/pm_distPower/getPowerLogCDF/main.F90
    !>  \compilef{getPowerLogCDF}
    !>  \output{getPowerLogCDF}
    !>  \include{lineno} example/pm_distPower/getPowerLogCDF/main.out.F90
    !>  \postproc{getPowerLogCDF}
    !>  \include{lineno} example/pm_distPower/getPowerLogCDF/main.py
    !>  \vis{getPowerLogCDF}
    !>  \image html pm_distPower/getPowerLogCDF/getPowerLogCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowerLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowerLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogCDFALD_RK5(logx, alpha, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogCDFALD_RK4(logx, alpha, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogCDFALD_RK3(logx, alpha, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogCDFALD_RK2(logx, alpha, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogCDFALD_RK1(logx, alpha, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, alpha, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogCDFALL_RK5(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogCDFALL_RK4(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogCDFALL_RK3(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogCDFALL_RK2(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogCDFALL_RK1(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogCDFALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Cumulative Distribution Function (CDF) of the
    !>  (Truncated) Power distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.
    !>
    !>  \param[out] logCDF      :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the CDF of the distribution.
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments, of type `real` of kind \RKALL,
    !>                              containing the natural logarithm of the value at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$)
    !>  \param[in]  logCDFNF    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Power distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getPowerLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPower::getPowerLogCDFNF).<br>
    !>
    !>  \interface{setPowerLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: setPowerLogCDF
    !>
    !>      call setPowerLogCDF(logCDF, logx, alpha, logCDFNF)
    !>      call setPowerLogCDF(logCDF, logx, alpha, logMinX, logCDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha` must hold for the corresponding input arguments.<br>
    !>  The condition `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Be warned that the procedures under this generic interface have no mechanism for checking for the consistency of the specified input arguments with each other.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowerLogCDF](@ref pm_distPower::setPowerLogCDF)<br>
    !>
    !>  \example{setPowerLogCDF}
    !>  \include{lineno} example/pm_distPower/setPowerLogCDF/main.F90
    !>  \compilef{setPowerLogCDF}
    !>  \output{setPowerLogCDF}
    !>  \include{lineno} example/pm_distPower/setPowerLogCDF/main.out.F90
    !>  \postproc{setPowerLogCDF}
    !>  \include{lineno} example/pm_distPower/setPowerLogCDF/main.py
    !>  \vis{setPowerLogCDF}
    !>  \image html pm_distPower/setPowerLogCDF/setPowerLogCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setPowerLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowerLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogCDFALD_RK5(logCDF, logx, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogCDFALD_RK4(logCDF, logx, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogCDFALD_RK3(logCDF, logx, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogCDFALD_RK2(logCDF, logx, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogCDFALD_RK1(logCDF, logx, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogCDFALL_RK5(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogCDFALL_RK4(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogCDFALL_RK3(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogCDFALL_RK2(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogCDFALL_RK1(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogCDFALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of the natural logarithm(s) of quantile corresponding to
    !>  the specified CDF of <b>(Truncated) Power distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the quantile of the (Truncated) Power distribution.<br>
    !>
    !>  \param[in]  logCDF      :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It can be present <b>if and only if</b> `logMaxX` is also present.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logx`                  :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the natural logarithm of the quantile corresponding to the input `logCDF`.<br>
    !>                              By definition, the condition `logMinX <= logx <= logMaxX` holds.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: getPowerLogQuan
    !>
    !>      logx = getPowerLogQuan(logCDF, alpha, logMaxX) ! Power distribution.
    !>      logx = getPowerLogQuan(logCDF, alpha, logMinX, logMaxX) ! Truncated Power distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logCDF < 0` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowerLogQuan](@ref pm_distPower::setPowerLogQuan)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distPower/getPowerLogQuan/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distPower/getPowerLogQuan/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distPower/getPowerLogQuan/main.py
    !>  \vis
    !>  \image html pm_distPower/getPowerLogQuan/getPowerLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowerLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogQuanALD_RK5(logCDF, alpha, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logCDF, alpha, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogQuanALD_RK4(logCDF, alpha, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logCDF, alpha, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogQuanALD_RK3(logCDF, alpha, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logCDF, alpha, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogQuanALD_RK2(logCDF, alpha, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logCDF, alpha, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogQuanALD_RK1(logCDF, alpha, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logCDF, alpha, logMaxX
        real(RKC)                               :: logx
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowerLogQuanALL_RK5(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowerLogQuanALL_RK4(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowerLogQuanALL_RK3(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowerLogQuanALL_RK2(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowerLogQuanALL_RK1(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogQuanALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the natural logarithm(s) of quantile corresponding to
    !>  the specified CDF of <b>(Truncated) Power distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the quantile of the (Truncated) Power distribution.<br>
    !>
    !>  \param[out] logx        :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the natural logarithm of the quantile corresponding to the input `logCDF`.<br>
    !>                              By definition, the condition `logMinX <= logx <= logMaxX` holds.<br>
    !>  \param[in]  logCDF      :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = `0`. It can be present <b>if and only if</b> `logCDFNF` is also present.)
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Power distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getPowerLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPower::getPowerLogCDFNF).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: setPowerLogQuan
    !>
    !>      call setPowerLogQuan(logx, logCDF, alpha, logCDFNF) ! Power distributed.
    !>      call setPowerLogQuan(logx, logCDF, alpha, logMinX, logCDFNF) ! Truncated Power distributed.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logCDF < 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowerLogQuan](@ref pm_distPower::getPowerLogQuan)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distPower/setPowerLogQuan/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distPower/setPowerLogQuan/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distPower/setPowerLogQuan/main.py
    !>  \vis
    !>  \image html pm_distPower/setPowerLogQuan/setPowerLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowerLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALD_RK5(logx, logCDF, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALD_RK4(logx, logCDF, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALD_RK3(logx, logCDF, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALD_RK2(logx, logCDF, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALD_RK1(logx, logCDF, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALL_RK5(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALL_RK4(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALL_RK3(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALL_RK2(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogQuanLLALL_RK1(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogQuanLLALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of the natural logarithm(s) of random value(s) from the
    !>  <b>(Truncated) Power distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more information on the (Truncated) Power distribution.<br>
    !>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = `0`. It can be present <b>if and only if</b> `logMaxX` is also present.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logRand`               :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the random value(s) from the specified distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: getPowerLogRand
    !>
    !>      logRand = getPowerLogRand(alpha, logMinX) ! Power distribution.
    !>      logRand = getPowerLogRand(alpha, logMinX, logMaxX) ! Truncated Power distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha > 0` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowerLogRand](@ref pm_distPower::setPowerLogRand)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distPower/getPowerLogRand/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distPower/getPowerLogRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distPower/getPowerLogRand/main.py
    !>  \vis
    !>  \image html pm_distPower/getPowerLogRand/getPowerLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowerLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getPowerLogRandALD_RK5(alpha, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getPowerLogRandALD_RK4(alpha, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getPowerLogRandALD_RK3(alpha, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getPowerLogRandALD_RK2(alpha, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getPowerLogRandALD_RK1(alpha, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getPowerLogRandALL_RK5(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getPowerLogRandALL_RK4(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getPowerLogRandALL_RK3(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getPowerLogRandALL_RK2(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getPowerLogRandALL_RK1(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowerLogRandALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the natural logarithm(s) of random value(s) from the
    !>  <b>(Truncated) Power distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  See the documentation of [pm_distPower](@ref pm_distPower) for more
    !>  information on generating random numbers from the (Truncated) Power distribution.
    !>
    !>  \param[out] logRand     :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the random value(s) from the specified distribution.<br>
    !>  \param[in]  negExpRand  :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing a random value from the standard Negative Exponential distribution (\f$\mu = 0, \sigma = 1.\f$).<br>
    !>                              This argument can be readily obtained by calling [getNegExpRand(sigma = 1.)](@ref pm_distNegExp::getNegExpRand).<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It can be present <b>if and only if</b> `logCDFNF` is also present.)
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Power distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getPowerLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPower::getPowerLogCDFNF).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distPower, only: setPowerLogRand
    !>
    !>      call setPowerLogRand(logRand, negExpRand, alpha, logCDFNF) ! Power distributed.
    !>      call setPowerLogRand(logRand, negExpRand, alpha, logMinX, logCDFNF) ! Truncated Power distributed.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `negExpRand <= 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowerLogRand](@ref pm_distPower::getPowerLogRand)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distPower/setPowerLogRand/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distPower/setPowerLogRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distPower/setPowerLogRand/main.py
    !>  \vis
    !>  \image html pm_distPower/setPowerLogRand/setPowerLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPower](@ref test_pm_distPower)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This interface can be extended to support vector-like `logRand` arguments other than the `elemental` approach.<br>
    !>  Such an extension would be sensible only if the new interface improves the performance against the `elemental` approach.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowerLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALD_RK5(logRand, negExpRand, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALD_RK4(logRand, negExpRand, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALD_RK3(logRand, negExpRand, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALD_RK2(logRand, negExpRand, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALD_RK1(logRand, negExpRand, alpha, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALL_RK5(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALL_RK4(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALL_RK3(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALL_RK2(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowerLogRandLNALL_RK1(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowerLogRandLNALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distPower
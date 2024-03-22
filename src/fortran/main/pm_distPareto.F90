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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>(Truncated) Pareto distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>(Truncated) Pareto distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function **(ICDF)** or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the <b>(Truncated) Pareto distribution</b> over a strictly-positive support \f$x \in [x_\mathrm{min}, x_\mathrm{max}]\f$
!>  is defined with the three <b>(shape, scale, scale)</b> parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ as,<br>
!>  \f{equation}{
!>      \large
!>      \pi(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~ x^{\alpha - 1} ~,
!>  \f}
!>  where \f$\mathbf{-\infty < \alpha < 0}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\eta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogPDFNF) of the PDF,<br>
!>  \f{equation}{
!>      \large
!>      \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~,
!>  \f}
!>
!>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the <b>Truncated Pareto distribution</b> simplifies to the **Pareto Distribution** with PDF,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{max} \rightarrow +\infty} \pi(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{-\alpha}{x_\mathrm{min}^\alpha} x^{\alpha - 1} ~,
!>  \f}
!>
!>  The equation for \f$\eta(\cdot)\f$ for the Pareto distribution simplifies to,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{max} \rightarrow +\infty} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = -\alpha x_\mathrm{min}^{-\alpha} ~,~ 0 < \alpha < +\infty.
!>  \f}
!>
!>  The corresponding **CDF** of the (Truncated) Pareto distribution is given by,<br>
!>  \f{equation}{
!>      \large
!>      \mathrm{CDF}(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) \bigg[1 - \bigg(\frac{x}{x_\mathrm{min}}\bigg)^\alpha\bigg] ~,
!>  \f}
!>
!>  where \f$\mathbf{0 < \alpha < +\infty}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogCDFNF) of the CDF,<br>
!>
!>  \f{eqnarray}{
!>      \large
!>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})
!>      &=& \frac{x_\mathrm{min}^\alpha}{x_\mathrm{min}^\alpha - x_\mathrm{max}^\alpha} ~, \\
!>      &=& -\frac{x_\mathrm{min}^\alpha}{\alpha} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~,
!>  \f}
!>
!>  where \f$\eta(\cdot)\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogPDFNF) of the PDF.<br>
!>
!>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the <b>Truncated Pareto distribution</b> simplifies to the **Pareto Distribution** with CDF,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{max} \rightarrow +\infty} \mathrm{CDF}(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \bigg[1 - \bigg(\frac{x}{x_\mathrm{min}}\bigg)^\alpha\bigg] ~,
!>  \f}
!>  with,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{x_\mathrm{max} \rightarrow +\infty} \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = +1 ~.
!>  \f}
!>
!>  The corresponding **Inverse CDF** of the (Truncated) Pareto distribution is given by,<br>
!>  \f{equation}{
!>      \large
!>      Q(\mathrm{CDF}(x); \alpha, x_\mathrm{min}, x_\mathrm{max}) \equiv x = x_\mathrm{min} \bigg(1 - \frac{\mathrm{CDF}(x)}{\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})}\bigg)^{\frac{1}{\alpha}} ~,
!>  \f}
!>  where \f$\mathbf{-\infty < \alpha < 0}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq Q(\mathrm{CDF}(x); \alpha, x_\mathrm{min}, x_\mathrm{max}) \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogCDFNF) of the CDF,<br>
!>  \f{equation}{
!>      \large
!>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{x_\mathrm{min}^\alpha}{x_\mathrm{min}^\alpha - x_\mathrm{max}^\alpha} ~.
!>  \f}
!>
!>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the **Truncated Pareto distribution** simplifies to the **Pareto Distribution**, with,<br>
!>  \f{equation}{
!>      \large
!>      Q(\mathrm{CDF}(x); \alpha, x_\mathrm{min}, x_\mathrm{max}) \equiv x = x_\mathrm{min} \big(1 - \mathrm{CDF}(x)\big)^{\frac{1}{\alpha}} ~,
!>  \f}
!>
!>  **Random Number Generation**<br>
!>
!>  Assuming that \f$U \in [0, 1)\f$ is a uniformly-distributed random variate, the transformed random variable,<br>
!>  \f{equation}{
!>      \large
!>      x = x_\mathrm{min} \bigg(1 - \frac{U}{\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})}\bigg)^{\frac{1}{\alpha}} ~,
!>  \f}
!>  follows a **Truncated Pareto distribution** with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max}\f$)
!>  where \f$\mathbf{-\infty < \alpha < 0}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogCDFNF) of the CDF,<br>
!>  \f{equation}{
!>      \large
!>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{x_\mathrm{min}^\alpha}{x_\mathrm{min}^\alpha - x_\mathrm{max}^\alpha} ~,
!>  \f}
!>
!>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the **Truncated Pareto distribution** simplifies to the **Pareto Distribution**, with,<br>
!>  \f{equation}{
!>      \large
!>      x = x_\mathrm{min} \big(1 - U\big)^{\frac{1}{\alpha}} ~,
!>  \f}
!>
!>  \remark
!>  See [pm_distPower](@ref pm_distPower) for the case of \f$0 < \alpha < +\infty\f$.<br>
!>  See [pm_distPoweto](@ref pm_distPoweto) for the case of \f$\alpha = 0\f$.<br>
!>
!>  \see
!>  [pm_distPower](@ref pm_distPower)<br>
!>  [pm_distPareto](@ref pm_distPareto)<br>
!>  [pm_distPoweto](@ref pm_distPoweto)<br>
!>
!>  \test
!>  [test_pm_distPareto](@ref test_pm_distPareto)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distPareto

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distPareto"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the PDF
    !>  of the (Truncated) Pareto distribution for an input parameter set \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  The normalization factor \f$\eta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ of the PDF of the Pareto distribution is,
    !>
    !>  \f{equation}{
    !>
    !>      \large
    !>      \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = \frac{\alpha}{x_\mathrm{max}^\alpha - x_\mathrm{min}^\alpha} ~,~ 0 < \alpha < +\infty
    !>
    !>  \f}
    !>
    !>  When \f$x_\mathrm{max} \rightarrow +\infty\f$, the <b>Truncated Pareto distribution</b> simplifies to the **Pareto Distribution**.<br>
    !>  The equation for \f$\eta(\cdot)\f$ for the Pareto distribution simplifies to,<br>
    !>
    !>  \f{equation}{
    !>      \large
    !>      \lim_{x_\mathrm{max} \rightarrow +\infty} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = -\alpha x_\mathrm{min}^{-\alpha} ~,~ 0 < \alpha < +\infty.
    !>  \f}
    !>
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factor of the PDF of the (Truncated) Pareto distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the properties of the (Truncated) Pareto distribution to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha   :   The input scalar or array of the same shape as other array-like arguments, of type `real` of kind \RKALL,
    !>                          containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>                          (**optional**, default = \f$+\infty\f$)
    !>
    !>  \return
    !>  `logPDFNF`          :   The output scalar or array of the same shape as any input array-like argument,
    !>                          of the same type and kind as the input argument `alpha`, containing the natural logarithm
    !>                          of the normalization factor of the PDF of (Truncated) Pareto distribution.<br>
    !>
    !>  \interface{getParetoLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: getParetoLogPDFNF
    !>
    !>      logPDFNF = getParetoLogPDFNF(alpha, logMinX) ! Pareto distribution.
    !>      logPDFNF = getParetoLogPDFNF(alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getParetoLogPDF](@ref pm_distPareto::getParetoLogPDF)<br>
    !>  [setParetoLogPDF](@ref pm_distPareto::setParetoLogPDF)<br>
    !>
    !>  \example{getParetoLogPDFNF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogPDFNF/main.F90
    !>  \compilef{getParetoLogPDFNF}
    !>  \output{getParetoLogPDFNF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogPDFNF/main.out.F90
    !>  \postproc{getParetoLogPDFNF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogPDFNF/main.py
    !>  \vis{getParetoLogPDFNF}
    !>  \image html pm_distPareto/getParetoLogPDFNF/getParetoLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getParetoLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getParetoLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogPDFNFALD_RK5(alpha, logMinX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogPDFNFALD_RK4(alpha, logMinX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogPDFNFALD_RK3(alpha, logMinX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogPDFNFALD_RK2(alpha, logMinX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogPDFNFALD_RK1(alpha, logMinX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogPDFNFALL_RK5(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogPDFNFALL_RK4(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogPDFNFALL_RK3(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogPDFNFALL_RK2(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogPDFNFALL_RK1(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFNFALL_RK1
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
    !>  (Truncated) Pareto distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `logx`, containing the natural logarithm of the PDF of the distribution.<br>
    !>
    !>  \interface{getParetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: getParetoLogPDF
    !>
    !>      logPDF = getParetoLogPDF(logx, alpha, logMinX) ! Pareto distribution.
    !>      logPDF = getParetoLogPDF(logx, alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  The conditions `logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setParetoLogPDF](@ref pm_distPareto::setParetoLogPDF)<br>
    !>
    !>  \example{getParetoLogPDF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogPDF/main.F90
    !>  \compilef{getParetoLogPDF}
    !>  \output{getParetoLogPDF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogPDF/main.out.F90
    !>  \postproc{getParetoLogPDF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogPDF/main.py
    !>  \vis{getParetoLogPDF}
    !>  \image html pm_distPareto/getParetoLogPDF/getParetoLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getParetoLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getParetoLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogPDFALD_RK5(logx, alpha, logMinX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogPDFALD_RK4(logx, alpha, logMinX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogPDFALD_RK3(logx, alpha, logMinX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogPDFALD_RK2(logx, alpha, logMinX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogPDFALD_RK1(logx, alpha, logMinX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogPDFALL_RK5(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogPDFALL_RK4(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogPDFALL_RK3(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogPDFALL_RK2(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogPDFALL_RK1(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogPDFALL_RK1
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
    !>  (Truncated) Pareto distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logPDFNF    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor of the PDF of the (Truncated) Pareto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getParetoLogPDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogPDFNF).<br>
    !>
    !>  \interface{setParetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: setParetoLogPDF
    !>
    !>      call setParetoLogPDF(logPDF, logx, alpha, logPDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
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
    !>  [setParetoLogPDF](@ref pm_distPareto::setParetoLogPDF)<br>
    !>
    !>  \example{setParetoLogPDF}
    !>  \include{lineno} example/pm_distPareto/setParetoLogPDF/main.F90
    !>  \compilef{setParetoLogPDF}
    !>  \output{setParetoLogPDF}
    !>  \include{lineno} example/pm_distPareto/setParetoLogPDF/main.out.F90
    !>  \postproc{setParetoLogPDF}
    !>  \include{lineno} example/pm_distPareto/setParetoLogPDF/main.py
    !>  \vis{setParetoLogPDF}
    !>  \image html pm_distPareto/setParetoLogPDF/setParetoLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setParetoLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setParetoLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogPDF_RK5(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogPDF_RK4(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogPDF_RK3(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogPDF_RK2(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: logx, alpha, logPDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogPDF_RK1(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogPDF_RK1
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
    !>  of the (Truncated) Pareto distribution for an input parameter set \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  The normalization factor \f$\zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ of the CDF of the Truncated Pareto distribution is,
    !>
    !>  \f{eqnarray}{
    !>      \large
    !>      \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max})
    !>      &=& \frac{x_\mathrm{min}^\alpha}{x_\mathrm{min}^\alpha - x_\mathrm{max}^\alpha} ~, \\
    !>      &=& -\frac{x_\mathrm{min}^\alpha}{\alpha} \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~,
    !>  \f}
    !>
    !>  where \f$\eta(\cdot)\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogPDFNF) of the PDF.<br>
    !>  The equation for \f$\zeta(\cdot)\f$ for the Pareto distribution simplifies to,<br>
    !>
    !>  \f{equation}{
    !>      \large
    !>      \lim_{x_\mathrm{max} \rightarrow +\infty} \zeta(\alpha, x_\mathrm{min}, x_\mathrm{max}) = +1 ~.
    !>  \f}
    !>
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factor of the CDF of the (Truncated) Pareto distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the properties of the (Truncated) Pareto distribution to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha   :   The input scalar or array of the same shape as other array-like arguments, of type `real` of kind \RKALL,
    !>                          containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `alpha`,
    !>                          containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>                          (**optional**, default = \f$+\infty\f$)
    !>
    !>  \return
    !>  `logCDFNF`          :   The output scalar or array of the same shape as any input array-like argument,
    !>                          of the same type and kind as the input argument `alpha`, containing the natural logarithm
    !>                          of the normalization factor of the CDF of (Truncated) Pareto distribution.<br>
    !>
    !>  \interface{getParetoLogCDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: getParetoLogCDFNF
    !>
    !>      logCDFNF = getParetoLogCDFNF(alpha, logMinX) ! = 0: Pareto distribution.
    !>      logCDFNF = getParetoLogCDFNF(alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The normalization factor in the case of the Pareto distribution is trivially `+1.`.
    !>
    !>  \see
    !>  [getParetoLogCDF](@ref pm_distPareto::getParetoLogCDF)<br>
    !>  [setParetoLogCDF](@ref pm_distPareto::setParetoLogCDF)<br>
    !>
    !>  \example{getParetoLogCDFNF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogCDFNF/main.F90
    !>  \compilef{getParetoLogCDFNF}
    !>  \output{getParetoLogCDFNF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogCDFNF/main.out.F90
    !>  \postproc{getParetoLogCDFNF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogCDFNF/main.py
    !>  \vis{getParetoLogCDFNF}
    !>  \image html pm_distPareto/getParetoLogCDFNF/getParetoLogCDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getParetoLogCDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getParetoLogCDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogCDFNFALD_RK5(alpha, logMinX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogCDFNFALD_RK4(alpha, logMinX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogCDFNFALD_RK3(alpha, logMinX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogCDFNFALD_RK2(alpha, logMinX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogCDFNFALD_RK1(alpha, logMinX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logCDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogCDFNFALL_RK5(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogCDFNFALL_RK4(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogCDFNFALL_RK3(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogCDFNFALL_RK2(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logCDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogCDFNFALL_RK1(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFNFALL_RK1
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
    !>  (Truncated) Pareto distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind the input argument `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$)
    !>
    !>  \return
    !>  `logCDF`                :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `logx`, containing the natural logarithm of the CDF of the distribution.<br>
    !>
    !>  \interface{getParetoLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: getParetoLogCDF
    !>
    !>      logCDF = getParetoLogCDF(logx, alpha, logMinX) ! Pareto distribution.
    !>      logCDF = getParetoLogCDF(logx, alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  The condition `logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setParetoLogCDF](@ref pm_distPareto::setParetoLogCDF)<br>
    !>
    !>  \example{getParetoLogCDF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogCDF/main.F90
    !>  \compilef{getParetoLogCDF}
    !>  \output{getParetoLogCDF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogCDF/main.out.F90
    !>  \postproc{getParetoLogCDF}
    !>  \include{lineno} example/pm_distPareto/getParetoLogCDF/main.py
    !>  \vis{getParetoLogCDF}
    !>  \image html pm_distPareto/getParetoLogCDF/getParetoLogCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getParetoLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getParetoLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogCDFALD_RK5(logx, alpha, logMinX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogCDFALD_RK4(logx, alpha, logMinX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogCDFALD_RK3(logx, alpha, logMinX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogCDFALD_RK2(logx, alpha, logMinX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogCDFALD_RK1(logx, alpha, logMinX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, alpha, logMinX
        real(RKC)                               :: logCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogCDFALL_RK5(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogCDFALL_RK4(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogCDFALL_RK3(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogCDFALL_RK2(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, alpha, logMinX, logMaxX
        real(RKC)                               :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogCDFALL_RK1(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogCDFALL_RK1
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
    !>  (Truncated) Pareto distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.
    !>
    !>  \param[out] logCDF      :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the CDF of the distribution.
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments, of type `real` of kind \RKALL,
    !>                              containing the natural logarithm of the value at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logCDFNF    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Pareto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getParetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogCDFNF).<br>
    !>                              (**optional**, default = [getParetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogCDFNF). If present, it implies a <b>(Truncated) Pareto distribution</b>.)
    !>
    !>  \interface{setParetoLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: setParetoLogCDF
    !>
    !>      call setParetoLogCDF(logCDF, logx, alpha, logMinX)
    !>      call setParetoLogCDF(logCDF, logx, alpha, logMinX, logCDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
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
    !>  [setParetoLogCDF](@ref pm_distPareto::setParetoLogCDF)<br>
    !>
    !>  \example{setParetoLogCDF}
    !>  \include{lineno} example/pm_distPareto/setParetoLogCDF/main.F90
    !>  \compilef{setParetoLogCDF}
    !>  \output{setParetoLogCDF}
    !>  \include{lineno} example/pm_distPareto/setParetoLogCDF/main.out.F90
    !>  \postproc{setParetoLogCDF}
    !>  \include{lineno} example/pm_distPareto/setParetoLogCDF/main.py
    !>  \vis{setParetoLogCDF}
    !>  \image html pm_distPareto/setParetoLogCDF/setParetoLogCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setParetoLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setParetoLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogCDFALD_RK5(logCDF, logx, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogCDFALD_RK4(logCDF, logx, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogCDFALD_RK3(logCDF, logx, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogCDFALD_RK2(logCDF, logx, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogCDFALD_RK1(logCDF, logx, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogCDFALL_RK5(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogCDFALL_RK4(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogCDFALL_RK3(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogCDFALL_RK2(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logCDF
        real(RKC)   , intent(in)                    :: logx, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogCDFALL_RK1(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogCDFALL_RK1
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
    !>  the specified CDF of <b>(Truncated) Pareto distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the quantile of the (Truncated) Pareto distribution.<br>
    !>
    !>  \param[in]  logCDF      :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$)
    !>
    !>  \return
    !>  `logx`                  :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the natural logarithm of the quantile corresponding to the input `logCDF`.<br>
    !>                              By definition, the condition `logMinX <= logx <= logMaxX` holds.<br>
    !>
    !>  \interface{getParetoLogQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: getParetoLogQuan
    !>
    !>      logx = getParetoLogQuan(logCDF, alpha, logMinX) ! Pareto distribution.
    !>      logx = getParetoLogQuan(logCDF, alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logCDF < 0` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setParetoLogQuan](@ref pm_distPareto::setParetoLogQuan)<br>
    !>
    !>  \example{getParetoLogQuan}
    !>  \include{lineno} example/pm_distPareto/getParetoLogQuan/main.F90
    !>  \compilef{getParetoLogQuan}
    !>  \output{getParetoLogQuan}
    !>  \include{lineno} example/pm_distPareto/getParetoLogQuan/main.out.F90
    !>  \postproc{getParetoLogQuan}
    !>  \include{lineno} example/pm_distPareto/getParetoLogQuan/main.py
    !>  \vis{getParetoLogQuan}
    !>  \image html pm_distPareto/getParetoLogQuan/getParetoLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getParetoLogQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getParetoLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogQuanALD_RK5(logCDF, alpha, logMinX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX
        real(RKC)                               :: logx
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogQuanALD_RK4(logCDF, alpha, logMinX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX
        real(RKC)                               :: logx
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogQuanALD_RK3(logCDF, alpha, logMinX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX
        real(RKC)                               :: logx
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogQuanALD_RK2(logCDF, alpha, logMinX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX
        real(RKC)                               :: logx
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogQuanALD_RK1(logCDF, alpha, logMinX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX
        real(RKC)                               :: logx
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getParetoLogQuanALL_RK5(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getParetoLogQuanALL_RK4(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getParetoLogQuanALL_RK3(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getParetoLogQuanALL_RK2(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logCDF, alpha, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getParetoLogQuanALL_RK1(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogQuanALL_RK1
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
    !>  the specified CDF of <b>(Truncated) Pareto distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the quantile of the (Truncated) Pareto distribution.<br>
    !>
    !>  \param[out] logx        :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the natural logarithm of the quantile corresponding to the input `logCDF`.<br>
    !>                              By definition, the condition `logMinX <= logx <= logMaxX` holds.<br>
    !>  \param[in]  logCDF      :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Pareto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getParetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogCDFNF).<br>
    !>                              (**optional**, default = [getParetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogCDFNF).
    !>                              If present, it implies a <b>(Truncated) Pareto distribution</b>.)
    !>
    !>  \interface{setParetoLogQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: setParetoLogQuan
    !>
    !>      call setParetoLogQuan(logx, logCDF, alpha, logMinX) ! Pareto distributed.
    !>      call setParetoLogQuan(logx, logCDF, alpha, logMinX, logCDFNF) ! Truncated Pareto distributed.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logCDF < 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getParetoLogQuan](@ref pm_distPareto::getParetoLogQuan)<br>
    !>
    !>  \example{setParetoLogQuan}
    !>  \include{lineno} example/pm_distPareto/setParetoLogQuan/main.F90
    !>  \compilef{setParetoLogQuan}
    !>  \output{setParetoLogQuan}
    !>  \include{lineno} example/pm_distPareto/setParetoLogQuan/main.out.F90
    !>  \postproc{setParetoLogQuan}
    !>  \include{lineno} example/pm_distPareto/setParetoLogQuan/main.py
    !>  \vis{setParetoLogQuan}
    !>  \image html pm_distPareto/setParetoLogQuan/setParetoLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \finmain{setParetoLogQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setParetoLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALD_RK5(logx, logCDF, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALD_RK4(logx, logCDF, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALD_RK3(logx, logCDF, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALD_RK2(logx, logCDF, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALD_RK1(logx, logCDF, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALL_RK5(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALL_RK4(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALL_RK3(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALL_RK2(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: logCDF, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogQuanLLALL_RK1(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogQuanLLALL_RK1
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
    !>  <b>(Truncated) Pareto distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more information on the (Truncated) Pareto distribution.<br>
    !>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$)
    !>
    !>  \return
    !>  `logRand`               :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the random value(s) from the specified distribution.<br>
    !>
    !>  \interface{getParetoLogRand}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: getParetoLogRand
    !>
    !>      logRand = getParetoLogRand(alpha, logMinX) ! Pareto distribution.
    !>      logRand = getParetoLogRand(alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setParetoLogRand](@ref pm_distPareto::setParetoLogRand)<br>
    !>
    !>  \example{getParetoLogRand}
    !>  \include{lineno} example/pm_distPareto/getParetoLogRand/main.F90
    !>  \compilef{getParetoLogRand}
    !>  \output{getParetoLogRand}
    !>  \include{lineno} example/pm_distPareto/getParetoLogRand/main.out.F90
    !>  \postproc{getParetoLogRand}
    !>  \include{lineno} example/pm_distPareto/getParetoLogRand/main.py
    !>  \vis{getParetoLogRand}
    !>  \image html pm_distPareto/getParetoLogRand/getParetoLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getParetoLogRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getParetoLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getParetoLogRandALD_RK5(alpha, logMinX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logRand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getParetoLogRandALD_RK4(alpha, logMinX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logRand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getParetoLogRandALD_RK3(alpha, logMinX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logRand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getParetoLogRandALD_RK2(alpha, logMinX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logRand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getParetoLogRandALD_RK1(alpha, logMinX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: alpha, logMinX
        real(RKC)                               :: logRand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getParetoLogRandALL_RK5(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getParetoLogRandALL_RK4(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getParetoLogRandALL_RK3(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getParetoLogRandALL_RK2(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: alpha, logMinX, logMaxX
        real(RKC)                               :: logRand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getParetoLogRandALL_RK1(alpha, logMinX, logMaxX) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParetoLogRandALL_RK1
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
    !>  <b>(Truncated) Pareto distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  See the documentation of [pm_distPareto](@ref pm_distPareto) for more
    !>  information on generating random numbers from the (Truncated) Pareto distribution.
    !>
    !>  \param[out] logRand     :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the random value(s) from the specified distribution.<br>
    !>  \param[in]  negExpRand  :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing a random value from the standard Negative Exponential distribution (\f$\mu = 0, \sigma = 1.\f$).<br>
    !>                              This argument can be readily obtained by calling [getNegExpRand(sigma = 1.)](@ref pm_distNegExp::getNegExpRand).<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the shape parameter (\f$\alpha\f$) of the distribution.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Pareto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getParetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogCDFNF).<br>
    !>                              (**optional**, default = [getParetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPareto::getParetoLogCDFNF).
    !>                              If present, it implies a <b>(Truncated) Pareto distribution</b>.)
    !>
    !>  \interface{setParetoLogRand}
    !>  \code{.F90}
    !>
    !>      use pm_distPareto, only: setParetoLogRand
    !>
    !>      call setParetoLogRand(logRand, negExpRand, alpha, logMinX) ! Pareto distributed.
    !>      call setParetoLogRand(logRand, negExpRand, alpha, logMinX, logCDFNF) ! Truncated Pareto distributed.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha < 0` must hold for the corresponding input arguments.<br>
    !>  The condition `negExpRand <= 0` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getParetoLogRand](@ref pm_distPareto::getParetoLogRand)<br>
    !>
    !>  \example{setParetoLogRand}
    !>  \include{lineno} example/pm_distPareto/setParetoLogRand/main.F90
    !>  \compilef{setParetoLogRand}
    !>  \output{setParetoLogRand}
    !>  \include{lineno} example/pm_distPareto/setParetoLogRand/main.out.F90
    !>  \postproc{setParetoLogRand}
    !>  \include{lineno} example/pm_distPareto/setParetoLogRand/main.py
    !>  \vis{setParetoLogRand}
    !>  \image html pm_distPareto/setParetoLogRand/setParetoLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPareto](@ref test_pm_distPareto)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This interface can be extended to support vector-like `logRand` arguments other than the `elemental` approach.<br>
    !>  Such an extension would be sensible only if the new interface improves the performance against the `elemental` approach.<br>
    !>
    !>  \finmain{setParetoLogRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setParetoLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALD_RK5(logRand, negExpRand, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALD_RK4(logRand, negExpRand, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALD_RK3(logRand, negExpRand, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALD_RK2(logRand, negExpRand, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALD_RK1(logRand, negExpRand, alpha, logMinX)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALL_RK5(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALL_RK4(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALL_RK3(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALL_RK2(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setParetoLogRandLNALL_RK1(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setParetoLogRandLNALL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha, logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distPareto
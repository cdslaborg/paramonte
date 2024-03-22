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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Beta distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Beta distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function **(ICDF)** or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the Beta distribution is defined with two shape parameters \f$(\alpha > 0, \beta > 0)\f$ as,
!>  \f{equation}{
!>      \large
!>      \pi(x | \alpha, \beta) =
!>      \begin{cases}
!>          \frac {x^{\alpha - 1} (1 - x)^{\beta - 1}} {\mathrm{B}(\alpha, \beta)} &,~ x \in (0, 1) ~, \\
!>          0 &,~ x \notin (0,1) ~.
!>      \end{cases}
!>  \f}
!>
!>  where,
!>  \f{equation}{
!>      \mathrm{B}(\alpha, \beta) = \frac {\Gamma(\alpha) \Gamma(\beta)} {\Gamma(\alpha + \beta)}
!>  \f}
!>
!>  is the Beta function whose natural logarithm is returned by [getLogBeta](@ref pm_mathBeta::getLogBeta).<br>
!>  Note that the common definitions of the Beta distribution consider the end points \f$ x = 0, x = 1\f$
!>  also as part of the support of the distribution.<br>
!>  But this inclusion is valid only if the condition \f$(\alpha \geq 1, \beta \geq 1)\f$ holds.<br>
!>  Furthermore, the end points correspond to a PDF value of either \f$0\f$ or \f$+\infty\f$,
!>  yielding a `NaN` when the natural logarithm of the PDF is returned in numerical algorithms.<br>
!>  As such, <b>the support of the Beta distribution in this module is always taken to be the open-interval</b> \f$x \in (0,1)\f$.<br>
!>
!>  The **CDF** of the Beta distribution is exactly given by the [regularized incomplete Beta function](@ref pm_mathBeta).<br>
!>  See the documentation of [pm_mathBeta](@ref pm_mathBeta) for details of the CDF of the Beta distribution.<br>
!>
!>  <b>Random Number Generation</b><br>
!>
!>  If \f$X\f$ and \f$Y\f$ are independent, with \f$X \sim \Gamma(\alpha, \theta)\f$ and \f$Y \sim \Gamma(\beta, \theta)\f$ then,
!>  \f{equation}{
!>      \frac{X}{X + Y} \sim \mathrm{B}(\alpha, \beta) ~,
!>  \f}
!>  where \f$(\alpha, \beta)\f$ are the parameters of the target [Beta distribution](@ref pm_distBeta) and,
!>  \f$\theta\f$ is an arbitrary parameter of the [Gamma distribution](@ref pm_distGamma).<br>
!>  So one algorithm for generating random Beta variates is to generate \f$\frac{X}{X + Y}\f$,
!>  where \f$X\f$ is a [Gamma variate](@ref pm_distGamma) with parameters \f$(\alpha, 1)\f$ and,
!>  \f$Y\f$ is an independent gamma variate with parameters \f$(\beta, 1)\f$.<br>
!>
!>  \note
!>  When \f$\alpha > 1\f$, the PDF of the Beta distribution is bound to zero at \f$x = 0\f$.<br>
!>  When \f$\beta  > 1\f$, the PDF of the Beta distribution is bound to zero at \f$x = 1\f$.<br>
!>  When \f$\alpha \leq 1\f$, the PDF of the Beta distribution is unbounded (\f$+\infty\f$) at \f$x = 0\f$.<br>
!>  When \f$\beta \leq 1\f$, the PDF of the Beta distribution is unbounded (\f$+\infty\f$) at \f$x = 1\f$.<br>
!>
!>  \note
!>  The **Cumulative Distribution Function (CDF)** of the Beta distribution for a triple \f$(x, \alpha, \beta)\f$
!>  is directly returned by calling [getBetaInc(x, alpha, beta)](@ref pm_mathBeta::getBetaInc).<br>
!>  Similarly, the **Cumulative Distribution Function (CDF)** of the Beta distribution for a triple \f$(x, \alpha, \beta)\f$
!>  is directly returned by calling [getBetaInv(x, alpha, beta)](@ref pm_mathBeta::getBetaInv).<br>
!>
!>  \see
!>  [pm_mathBeta](@ref pm_mathBeta)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distGamma](@ref pm_distGamma)<br>
!>  [pm_mathGamma](@ref pm_mathGamma)<br>
!>
!>  \test
!>  [test_pm_distBeta](@ref test_pm_distBeta)
!>
!>  \todo
!>  \pvhigh
!>  The quantile function must be implemented here.
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distBeta

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distBeta"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Probability Density Function (PDF) of the
    !>  Beta distribution for an input `x` within the support of the distribution \f$x \in [0,1]\f$.
    !>
    !>  \brief
    !>  See the documentation of [pm_distBeta](@ref pm_distBeta) for more information on the Beta distribution.<br>
    !>
    !>  Contrary to the higher-performance procedure interfaces [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF) and
    !>  [setBetaLogPDF](@ref pm_distBeta::setBetaLogPDF) that only accept \f$x \in (0,1)\f$,
    !>  the procedures under this generic interface also accept the boundary values of the support: \f$x = 0\f$, and \f$x = 1\f$.<br>
    !>
    !>  This interface therefore has two primary use cases:
    !>  <ol>
    !>      <li>    When the PDF of the Beta distribution (as opposed to the natural logarithm of PDF) is needed.
    !>      <li>    When the coverage of the full closed interval \f$[0,1]\f$ for \f$x\f$ is desired.
    !>  </ol>
    !>  Note, however, the output of the procedures of this generic interface are
    !>  **prone to overflow** at \f$(x \sim 0, \alpha < 1)\f$ or \f$(x \sim 1, \beta < 1)\f$.<br>
    !>  The occurrence of the overflow depends on the specified values for \f$(x, \alpha, \beta)\f$
    !>  and the requested `real` precision (kind) for the output PDF value.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the second shape parameter of the distribution.<br>
    !>
    !>  \return
    !>  `pdf`                   :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `x`, containing the PDF of the distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distBeta, only: getBetaPDF
    !>
    !>      pdf = getBetaPDF(x, alpha, beta)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition \f$x \in (0,1)\f$, \f$\alpha \geq 0\f$, and \f$\beta \geq 0\f$ must hold for the triple input values `(x, alpha, beta)`.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The procedures of this generic interface are prone to overflow at extreme values of `x, alpha, beta`.<br>
    !>  Use [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF) or [setBetaLogPDF](@ref pm_distBeta::setBetaLogPDF)
    !>  which more tolerant of extreme input values.
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>  [setBetaLogPDF](@ref pm_distBeta::setBetaLogPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distBeta/getBetaPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distBeta/getBetaPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distBeta/getBetaPDF/main.py
    !>  \vis
    !>  \image html pm_distBeta/getBetaPDF/getBetaPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBeta](@ref test_pm_distBeta)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBetaPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBetaPDF_RK5(x, alpha, beta) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: pdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBetaPDF_RK4(x, alpha, beta) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: pdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBetaPDF_RK3(x, alpha, beta) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: pdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBetaPDF_RK2(x, alpha, beta) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: pdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBetaPDF_RK1(x, alpha, beta) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: pdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the
    !>  Beta distribution for an input `x` within the support of the distribution \f$x \in (0,1)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distBeta](@ref pm_distBeta) for more information on the Beta distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the second shape parameter of the distribution.<br>
    !>  \param[in]  logBeta     :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the natural logarithm of the Beta function \f$\mathrm{B}(\alpha, \beta)\f$.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, \beta)\f$ parameter
    !>                              will significantly improve the runtime performance.<br>
    !>                              (**optional**, default = [getLogBeta(alpha, beta)](@ref pm_mathBeta::getLogBeta))
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distBeta, only: getBetaLogPDF
    !>
    !>      logPDF = getBetaLogPDF(x, alpha, beta)
    !>      logPDF = getBetaLogPDF(x, alpha, beta, logBeta)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition \f$x \in (0,1)\f$, \f$\alpha > 0\f$, and \f$\beta > 0\f$ must hold for the triple input values `(x, alpha, beta)`.<br>
    !>  Also, the value of the input argument `logBeta`, if present, must be consistent with input values for `alpha` and `beta`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setBetaLogPDF](@ref pm_distBeta::setBetaLogPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distBeta/getBetaLogPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distBeta/getBetaLogPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distBeta/getBetaLogPDF/main.py
    !>  \vis
    !>  \image html pm_distBeta/getBetaLogPDF/getBetaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBeta](@ref test_pm_distBeta)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBetaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBetaLogPDFD_RK5(x, alpha, beta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBetaLogPDFD_RK4(x, alpha, beta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBetaLogPDFD_RK3(x, alpha, beta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBetaLogPDFD_RK2(x, alpha, beta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBetaLogPDFD_RK1(x, alpha, beta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: x, alpha, beta
        real(RKC)                   :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBetaLogPDFL_RK5(x, alpha, beta, logBeta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: x, alpha, beta, logBeta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBetaLogPDFL_RK4(x, alpha, beta, logBeta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: x, alpha, beta, logBeta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBetaLogPDFL_RK3(x, alpha, beta, logBeta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: x, alpha, beta, logBeta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBetaLogPDFL_RK2(x, alpha, beta, logBeta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: x, alpha, beta, logBeta
        real(RKC)                   :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBetaLogPDFL_RK1(x, alpha, beta, logBeta) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaLogPDFL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: x, alpha, beta, logBeta
        real(RKC)                   :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the
    !>  Beta distribution for an input `x` within the support of the distribution \f$x \in (0,1)\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distBeta](@ref pm_distBeta) for more information on the Beta distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the values at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the second shape parameter of the distribution.<br>
    !>  \param[in]  logBeta     :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `x`, containing the natural logarithm of the Beta function \f$\mathrm{B}(\alpha, \beta)\f$.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, \beta)\f$ parameter
    !>                              will significantly improve the runtime performance.<br>
    !>                              (**optional**, default = [getLogBeta(alpha, beta)](@ref pm_mathBeta::getLogBeta))
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distBeta, only: setBetaLogPDF
    !>
    !>      call setBetaLogPDF(logPDF, x, alpha, beta)
    !>      call setBetaLogPDF(logPDF, x, alpha, beta, logBeta)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition \f$x \in (0,1)\f$, \f$\alpha > 0\f$, and \f$\beta > 0\f$ must hold for the triple input values `(x, alpha, beta)`.<br>
    !>  Also, the value of the input argument `logBeta`, if present, must be consistent with input values for `alpha` and `beta`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setBetaLogPDF](@ref pm_distBeta::setBetaLogPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distBeta/setBetaLogPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distBeta/setBetaLogPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distBeta/setBetaLogPDF/main.py
    !>  \vis
    !>  \image html pm_distBeta/setBetaLogPDF/setBetaLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBeta](@ref test_pm_distBeta)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBetaLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBetaLogPDFD_RK5(logPDF, x, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBetaLogPDFD_RK4(logPDF, x, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBetaLogPDFD_RK3(logPDF, x, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBetaLogPDFD_RK2(logPDF, x, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBetaLogPDFD_RK1(logPDF, x, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBetaLogPDFL_RK5(logPDF, x, alpha, beta, logBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta, logBeta
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBetaLogPDFL_RK4(logPDF, x, alpha, beta, logBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta, logBeta
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBetaLogPDFL_RK3(logPDF, x, alpha, beta, logBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta, logBeta
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBetaLogPDFL_RK2(logPDF, x, alpha, beta, logBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta, logBeta
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBetaLogPDFL_RK1(logPDF, x, alpha, beta, logBeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaLogPDFL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, alpha, beta, logBeta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the CDF of the Beta distribution for the given parameters \f$(\alpha, \beta)\f$
    !>  as defined in the details section of [pm_distBeta](@ref pm_distBeta).
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
    !>                                                      Therefore, the user is expected to be aware of and apply the necessary correction (addition by `1`) before using the output value.<br>
    !>                                          </ol>
    !>                              </ol>
    !>                              Specifying `signed = .true.` can lead to considerably more accurate calculations (by orders of magnitudes) for values of `x` and `betaInc` that are near `1`.<br>
    !>                              The loss of precision near `1` occurs because of inadequate resolution of `real` number representations in digital computers near `1`
    !>                              which is orders of magnitude worse than the precision near `0`.<br>
    !>                              (**optional**, default = `.false.`, following the principle of least surprise.)
    !>
    !>  \return
    !>  `cdf`                   :   The output object of the same type, kind, rank, and shape as other input array-like arguments
    !>                              containing the CDF of the Beta distribution.
    !>
    !>  \interface{getBetaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distBeta, only: getBetaCDF
    !>
    !>      cdf = getBetaCDF(x, alpha, beta, signed = signed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings associated with [getBetaInc](@ref pm_mathBeta::getBetaInc)
    !>  also apply to the procedures of this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [getBetaCDF](@ref pm_distBeta::getBetaCDF)<br>
    !>  [setBetaRand](@ref pm_distBeta::setBetaRand)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{getBetaCDF}
    !>  \include{lineno} example/pm_distBeta/getBetaCDF/main.F90
    !>  \compilef{getBetaCDF}
    !>  \output{getBetaCDF}
    !>  \include{lineno} example/pm_distBeta/getBetaCDF/main.out.F90
    !>  \postproc{getBetaCDF}
    !>  \include{lineno} example/pm_distBeta/getBetaCDF/main.py
    !>  \vis{getBetaCDF}
    !>  \image html pm_distBeta/getBetaCDF/getBetaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBeta](@ref test_pm_distBeta)
    !>
    !>  \finmain{getBetaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBetaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getBetaCDF_RK5(x, alpha, beta, signed) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK) , intent(in), optional  :: signed
        real(RKC)   , intent(in)            :: x, alpha, beta
        real(RKC)                           :: cdf
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getBetaCDF_RK4(x, alpha, beta, signed) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK) , intent(in), optional  :: signed
        real(RKC)   , intent(in)            :: x, alpha, beta
        real(RKC)                           :: cdf
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getBetaCDF_RK3(x, alpha, beta, signed) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK) , intent(in), optional  :: signed
        real(RKC)   , intent(in)            :: x, alpha, beta
        real(RKC)                           :: cdf
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getBetaCDF_RK2(x, alpha, beta, signed) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK) , intent(in), optional  :: signed
        real(RKC)   , intent(in)            :: x, alpha, beta
        real(RKC)                           :: cdf
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getBetaCDF_RK1(x, alpha, beta, signed) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK) , intent(in), optional  :: signed
        real(RKC)   , intent(in)            :: x, alpha, beta
        real(RKC)                           :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the CDF of the Beta distribution for the given parameters \f$(\alpha, \beta)\f$
    !>  as defined in the details section of [pm_distBeta](@ref pm_distBeta).
    !>
    !>  \param[in]  cdf         :   The output scalar (or array of the same shape as other array-like arguments)
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
    !>                                                      Therefore, the user is expected to be aware of and apply the necessary correction (addition by `1`) before using the output value.<br>
    !>                                          </ol>
    !>                              </ol>
    !>                              Specifying `signed = .true.` can lead to considerably more accurate calculations (by orders of magnitudes) for values of `x` and `betaInc` that are near `1`.<br>
    !>                              The loss of precision near `1` occurs because of inadequate resolution of `real` number representations in digital computers near `1`
    !>                              which is orders of magnitude worse than the precision near `0`.<br>
    !>  \param[out] info        :   The output scalar (or array of the same shape as other array-like arguments) of type `integer` of default kind \IK.<br>
    !>                              It is `0` <b>if and only if</b> the computation converges, otherwise,
    !>                              it is set to a non-zero integer whose value depends on the root cause of the lack of convergence.<br>
    !>                              Convergence typically fails if `alpha` or `beta` are too large.<br>
    !>                              For more information, see the corresponding argument of [setBetaInc](@ref pm_mathBeta::setBetaInc).<br>
    !>
    !>  \return
    !>
    !>  \interface{setBetaCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distBeta, only: setBetaCDF
    !>
    !>      call setBetaCDF(cdf, x, alpha, beta, logFuncBeta, signed, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings associated with [setBetaInc](@ref pm_mathBeta::setBetaInc) also apply to the procedures of this generic interface.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogBeta](@ref pm_mathBeta::getLogBeta)<br>
    !>  [setBetaCDF](@ref pm_distBeta::setBetaCDF)<br>
    !>  [setBetaRand](@ref pm_distBeta::setBetaRand)<br>
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>
    !>  \example{setBetaCDF}
    !>  \include{lineno} example/pm_distBeta/setBetaCDF/main.F90
    !>  \compilef{setBetaCDF}
    !>  \output{setBetaCDF}
    !>  \include{lineno} example/pm_distBeta/setBetaCDF/main.out.F90
    !>  \postproc{setBetaCDF}
    !>  \include{lineno} example/pm_distBeta/setBetaCDF/main.py
    !>  \vis{setBetaCDF}
    !>  \image html pm_distBeta/setBetaCDF/setBetaCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBeta](@ref test_pm_distBeta)
    !>
    !>  \finmain{setBetaCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBetaCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBetaCDF_RK5(cdf, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)   :: cdf
        real(RKC)   , intent(in)    :: x, alpha, beta, logFuncBeta
        logical(LK) , intent(in)    :: signed
        integer(IK) , intent(out)   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBetaCDF_RK4(cdf, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)   :: cdf
        real(RKC)   , intent(in)    :: x, alpha, beta, logFuncBeta
        logical(LK) , intent(in)    :: signed
        integer(IK) , intent(out)   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBetaCDF_RK3(cdf, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)   :: cdf
        real(RKC)   , intent(in)    :: x, alpha, beta, logFuncBeta
        logical(LK) , intent(in)    :: signed
        integer(IK) , intent(out)   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBetaCDF_RK2(cdf, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)   :: cdf
        real(RKC)   , intent(in)    :: x, alpha, beta, logFuncBeta
        logical(LK) , intent(in)    :: signed
        integer(IK) , intent(out)   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBetaCDF_RK1(cdf, x, alpha, beta, logFuncBeta, signed, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)   :: cdf
        real(RKC)   , intent(in)    :: x, alpha, beta, logFuncBeta
        logical(LK) , intent(in)    :: signed
        integer(IK) , intent(out)   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar or array of arbitrary rank of Beta-distributed random values in range \f$[0, 1]\f$
    !>  (or \f$(0, 1)\f$, depending on the specific parameter values) with the specified two shape parameters
    !>  \f$(\alpha, \beta)\f$ of the Beta distribution corresponding to the procedure arguments `(alpha, beta)`.
    !>
    !>  \details
    !>  See the documentation of [pm_distBeta](@ref pm_distBeta) for more
    !>  information on the Probability Density Function (PDF) of the Beta distribution and RNG.
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type), implying the use of the intrinsic Fortran URNG.)
    !>  \param[out]     rand    :   The output scalar or
    !>                              <ol>
    !>                                  <li>    array of rank `1`, or<br>
    !>                                  <li>    array of arbitrary rank if the `rng` argument is missing or set to [rngf_type](@ref pm_distUnif::rngf_type), or<br>
    !>                              </ol>
    !>                              of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL.
    !>                              </ol>
    !>                              On output, it contains Beta-distributed random value(s).<br>
    !>  \param[in]      alpha   :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `rand`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]      beta    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `rand`, containing the second shape parameter of the distribution.<br>
    !>
    !>  \interface{setBetaRand}
    !>  \code{.F90}
    !>
    !>      use pm_distBeta, only: setBetaRand
    !>
    !>      call setBetaRand(rand, alpha, beta)
    !>      call setBetaRand(rand(..), alpha, beta)
    !>      call setBetaRand(rng, rand, alpha, beta)
    !>      call setBetaRand(rng, rand(:), alpha, beta)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions \f$0. < \alpha\f$ and \f$0. < \beta\f$ must hold for the input arguments `(alpha, beta)`.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \recursive
    !>
    !>  \note
    !>  For repeated Beta RNG with fixed `alpha`, it is best to pass a vector of `rand` to be filled
    !>  with random numbers rather than calling the procedures with scalar `rand` argument repeatedly.<br>
    !>  In addition to avoiding procedure call overhead, vectorized RGN in this particular case also avoids
    !>  an unnecessary division and square-root operation.<br>
    !>
    !>  \see
    !>  [getBetaLogPDF](@ref pm_distBeta::getBetaLogPDF)<br>
    !>  [setBetaLogPDF](@ref pm_distBeta::setBetaLogPDF)<br>
    !>  [getBetaCDF](@ref pm_distBeta::getBetaCDF)<br>
    !>  [setBetaCDF](@ref pm_distBeta::setBetaCDF)<br>
    !>
    !>  \example{setBetaRand}
    !>  \include{lineno} example/pm_distBeta/setBetaRand/main.F90
    !>  \compilef{setBetaRand}
    !>  \output{setBetaRand}
    !>  \include{lineno} example/pm_distBeta/setBetaRand/main.out.F90
    !>  \postproc{setBetaRand}
    !>  \include{lineno} example/pm_distBeta/setBetaRand/main.py
    !>  \vis{setBetaRand}
    !>  \image html pm_distBeta/setBetaRand/setBetaRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBeta](@ref test_pm_distBeta)
    !>
    !>  \finmain{setBetaRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBetaRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBetaRandRNGD_D0_RK5(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBetaRandRNGD_D0_RK4(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBetaRandRNGD_D0_RK3(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBetaRandRNGD_D0_RK2(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBetaRandRNGD_D0_RK1(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBetaRandRNGF_D0_RK5(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBetaRandRNGF_D0_RK4(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBetaRandRNGF_D0_RK3(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBetaRandRNGF_D0_RK2(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBetaRandRNGF_D0_RK1(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setBetaRandRNGX_D0_RK5(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setBetaRandRNGX_D0_RK4(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setBetaRandRNGX_D0_RK3(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setBetaRandRNGX_D0_RK2(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setBetaRandRNGX_D0_RK1(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setBetaRandRNGD_D1_RK5(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setBetaRandRNGD_D1_RK4(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setBetaRandRNGD_D1_RK3(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setBetaRandRNGD_D1_RK2(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setBetaRandRNGD_D1_RK1(rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGD_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setBetaRandRNGF_D1_RK5(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setBetaRandRNGF_D1_RK4(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setBetaRandRNGF_D1_RK3(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setBetaRandRNGF_D1_RK2(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setBetaRandRNGF_D1_RK1(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGF_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setBetaRandRNGX_D1_RK5(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setBetaRandRNGX_D1_RK4(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setBetaRandRNGX_D1_RK3(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setBetaRandRNGX_D1_RK2(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setBetaRandRNGX_D1_RK1(rng, rand, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBetaRandRNGX_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: alpha, beta
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distBeta
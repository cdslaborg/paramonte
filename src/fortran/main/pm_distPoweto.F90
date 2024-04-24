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
!>  This module contains classes and procedures for computing various statistical quantities
!>  related to the <b>(Truncated) Power/Pareto distribution</b> (hence the name **Poweto**).
!>
!>  \details
!>  Specifically, this module contains routines for computing the
!>  following quantities of the <b>(Truncated) Power/Pareto distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The (Truncated) Poweto distribution is a generalization of the (Truncated)
!>  [Pareto](@ref pm_distPareto) and [Power](@ref pm_distPower) distributions.<br>
!>
!>  The **PDF** of the <b>(Truncated) Poweto distribution</b> over a strictly-positive support \f$x \in [x_\mathrm{min}, x_\mathrm{max}]\f$
!>  is defined with the three <b>(shape, scale, scale)</b> parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ as,<br>
!>  \f{equation}{
!>      \large
!>      \pi(x | \alpha, x_\mathrm{min}, x_\mathrm{max}) = \eta(\alpha, x_\mathrm{min}, x_\mathrm{max}) ~ x^{\alpha - 1} ~,
!>  \f}
!>  where \f$\mathbf{-\infty < \alpha < +\infty}\f$ and \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ hold,
!>  and \f$\eta(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distPareto::getParetoLogPDFNF) of the PDF.<br>
!>
!>  When \f$\alpha = 0\f$, the <b>(Truncated) Poweto distribution</b> simplifies to the Power-law distribution with exponent \f$-1\f$.<br>
!>  When \f$0 < \alpha < +\infty\f$, the <b>(Truncated) Poweto distribution</b> simplifies to [(Truncated) Power distribution](@ref pm_distPower).<br>
!>  When \f$-\infty < \alpha < 0\f$, the <b>(Truncated) Poweto distribution</b> simplifies to [(Truncated) Pareto distribution](@ref pm_distPareto).<br>
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distPower](@ref pm_distPower)<br>
!>  [pm_distPareto](@ref pm_distPareto)<br>
!>  [pm_distPoweto](@ref pm_distPoweto)<br>
!>  [pm_distPiwiPoweto](@ref pm_distPiwiPoweto)<br>
!>
!>  \test
!>  [test_pm_distPoweto](@ref test_pm_distPoweto)
!>
!>  \todo
!>  Generic interfaces for computing the logarithm of CDF robustly (without numerical rounding) must be added in the future.
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distPoweto

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distPoweto"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the Probability Density Function (PDF)
    !>  of the (Truncated) Poweto distribution for the input parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for the definition of the normalization factor.<br>
    !>
    !>  The primary use of this interface is to compute the normalization
    !>  factors of the (Truncated) Poweto distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the (Truncated) Poweto PDF
    !>  to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of type `real` of kind \RKALL,
    !>                              containing the shape parameter of the distribution (i.e., the exponent plus one of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution,
    !>                              representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$. It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \return
    !>  `logPDFNF`              :   The output scalar of the same type, kind, and size as the input argument `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the power-law
    !>                              component of the (Truncated) Poweto distribution.<br>
    !>
    !>  \interface{getPowetoLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogPDFNF
    !>
    !>      logPDFNF = getPowetoLogPDFNF(alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowetoLogPDF](@ref pm_distPoweto::getPowetoLogPDF)<br>
    !>  [setPowetoLogPDF](@ref pm_distPoweto::setPowetoLogPDF)<br>
    !>
    !>  \example{getPowetoLogPDFNF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogPDFNF/main.F90
    !>  \compilef{getPowetoLogPDFNF}
    !>  \output{getPowetoLogPDFNF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogPDFNF/main.out.F90
    !>  \postproc{getPowetoLogPDFNF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogPDFNF/main.py
    !>  \vis
    !>  \image html pm_distPoweto/getPowetoLogPDFNF/getPowetoLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowetoLogPDFNF_D0_RK5(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowetoLogPDFNF_D0_RK4(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowetoLogPDFNF_D0_RK3(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowetoLogPDFNF_D0_RK2(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowetoLogPDFNF_D0_RK1(alpha, logMinX, logMaxX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of
    !>  the (Truncated) Poweto distribution for an input `logx` within the support of
    !>  the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar (or array of the same shape as other array-like arguments)
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the `x` value
    !>                              within the support of the distribution at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution,
    !>                              representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar of the same type and kind the input argument `logx`,
    !>                              containing the natural logarithm of the PDF of the distribution at the specified point within the support of the PDF.<br>
    !>
    !>  \interface{getPowetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogPDF
    !>
    !>      logPDF = getPowetoLogPDF(logx, alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX <= logx .and. logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPowetoLogPDF](@ref pm_distPoweto::setPowetoLogPDF)<br>
    !>
    !>  \example{getPowetoLogPDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogPDF/main.F90
    !>  \compilef{getPowetoLogPDF}
    !>  \output{getPowetoLogPDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogPDF/main.out.F90
    !>  \postproc{getPowetoLogPDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogPDF/main.py
    !>  \vis{getPowetoLogPDF}
    !>  \image html pm_distPoweto/getPowetoLogPDF/getPowetoLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowetoLogPDF_D0_RK5(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowetoLogPDF_D0_RK4(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowetoLogPDF_D0_RK3(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowetoLogPDF_D0_RK2(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowetoLogPDF_D0_RK1(logx, alpha, logMinX, logMaxX) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of
    !>  the (Truncated) Poweto distribution for an input `logx` within the support of
    !>  the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar of the same type and kind the input argument `logx`, containing the natural logarithm of the PDF of the distribution.<br>
    !>  \param[in]  logx        :   The input scalar (or array of the same shape as other array-like arguments)
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the `x` value
    !>                              within the support of the distribution at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logPDFNF    :   The input scalar (or array of the same shape as other array-like arguments) of the same type, kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor (\f$\eta\f$) of power-law function of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters significantly improves the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getPowetoLogPDFNF(alpha, logMinX, logMaxX)](@ref pm_distPoweto::getPowetoLogPDFNF).<br>
    !>
    !>  \interface{setPowetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: setPowetoLogPDF
    !>
    !>      call setPowetoLogPDF(logPDF, logx, alpha, logPDFNF)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowetoLogPDF](@ref pm_distPoweto::getPowetoLogPDF)<br>
    !>
    !>  \example{setPowetoLogPDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogPDF/main.F90
    !>  \compilef{setPowetoLogPDF}
    !>  \output{setPowetoLogPDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogPDF/main.out.F90
    !>  \postproc{setPowetoLogPDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogPDF/main.py
    !>  \vis{setPowetoLogPDF}
    !>  \image html pm_distPoweto/setPowetoLogPDF/setPowetoLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setPowetoLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowetoLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowetoLogPDF_D0_RK5(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                            :: logPDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowetoLogPDF_D0_RK4(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                            :: logPDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowetoLogPDF_D0_RK3(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                            :: logPDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowetoLogPDF_D0_RK2(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                            :: logPDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowetoLogPDF_D0_RK1(logPDF, logx, alpha, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                            :: logPDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factor of the Cumulative Distribution Function (CDF)
    !>  of the (Truncated) Poweto distribution for the input parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for the definition of the normalization factor.<br>
    !>
    !>  The primary use of this interface is to compute the normalization
    !>  factors of the (Truncated) Poweto distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the (Truncated) Poweto CDF
    !>  to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of type `real` of kind \RKALL,
    !>                              containing the shape parameter of the distribution (i.e., the exponent plus one of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution,
    !>                              representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \return
    !>  `logCDFNF`              :   The output scalar of the same type, kind, and size as the input argument `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the power-law
    !>                              component of the (Truncated) Poweto distribution.<br>
    !>
    !>  \interface{getPowetoLogCDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogCDFNF
    !>
    !>      logCDFNF = getPowetoLogCDFNF(alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>      logCDFNF = getPowetoLogCDFNF(alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowetoLogCDF](@ref pm_distPoweto::getPowetoLogCDF)<br>
    !>  [setPowetoLogCDF](@ref pm_distPoweto::setPowetoLogCDF)<br>
    !>
    !>  \example{getPowetoLogCDFNF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogCDFNF/main.F90
    !>  \compilef{getPowetoLogCDFNF}
    !>  \output{getPowetoLogCDFNF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogCDFNF/main.out.F90
    !>  \postproc{getPowetoLogCDFNF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogCDFNF/main.py
    !>  \vis
    !>  \image html pm_distPoweto/getPowetoLogCDFNF/getPowetoLogCDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogCDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogCDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowetoLogCDFNF_D0_RK5(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDFNF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logCDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowetoLogCDFNF_D0_RK4(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDFNF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logCDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowetoLogCDFNF_D0_RK3(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDFNF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logCDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowetoLogCDFNF_D0_RK2(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDFNF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logCDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowetoLogCDFNF_D0_RK1(alpha, logMinX, logMaxX) result(logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDFNF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX
        real(RKC)                                   :: logCDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Cumulative Distribution Function (CDF) of
    !>  the (Truncated) Poweto distribution for an input `logx` within the support of
    !>  the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar (or array of the same shape as other array-like arguments)
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the `x` value
    !>                              within the support of the distribution at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution,
    !>                              representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$. It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \return
    !>  `logCDF`                :   The output scalar of the same type and kind the input argument `logx`,
    !>                              containing the natural logarithm of the CDF of the distribution at the specified point within the distribution support.<br>
    !>
    !>  \interface{getPowetoLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogCDF
    !>
    !>      logCDF = getPowetoLogCDF(logx, alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX <= logx .and. logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPowetoLogCDF](@ref pm_distPoweto::setPowetoLogCDF)<br>
    !>
    !>  \example{getPowetoLogCDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogCDF/main.F90
    !>  \compilef{getPowetoLogCDF}
    !>  \output{getPowetoLogCDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogCDF/main.out.F90
    !>  \postproc{getPowetoLogCDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogCDF/main.py
    !>  \vis{getPowetoLogCDF}
    !>  \image html pm_distPoweto/getPowetoLogCDF/getPowetoLogCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowetoLogCDF_D0_RK5(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowetoLogCDF_D0_RK4(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowetoLogCDF_D0_RK3(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowetoLogCDF_D0_RK2(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowetoLogCDF_D0_RK1(logx, alpha, logMinX, logMaxX) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogCDF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Cumulative Distribution Function (CDF) of
    !>  the (Truncated) Poweto distribution for an input `logx` within the support of
    !>  the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[out] logCDF      :   The output scalar of the same type and kind the input argument `logx`, containing the natural logarithm of the CDF of the distribution.<br>
    !>  \param[in]  logx        :   The input scalar (or array of the same shape as other array-like arguments)
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the `x` value
    !>                              within the support of the distribution at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `logx`,
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same shape as other array-like arguments) of the same type, kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor (\f$\eta\f$) of power-law function of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters significantly improves the runtime performance.<br>
    !>                              (**optional**, default = [getPowetoLogCDFNF(alpha, logMinX)](@ref pm_distPareto::getParetoLogCDFNF). It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \interface{setPowetoLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: setPowetoLogCDF
    !>
    !>      call setPowetoLogCDF(logCDF, logx, alpha, logMinX = logMinX, logCDFNF = logCDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logCDFNF)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowetoLogCDF](@ref pm_distPoweto::getPowetoLogCDF)<br>
    !>
    !>  \example{setPowetoLogCDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogCDF/main.F90
    !>  \compilef{setPowetoLogCDF}
    !>  \output{setPowetoLogCDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogCDF/main.out.F90
    !>  \postproc{setPowetoLogCDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogCDF/main.py
    !>  \vis{setPowetoLogCDF}
    !>  \image html pm_distPoweto/setPowetoLogCDF/setPowetoLogCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setPowetoLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowetoLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowetoLogCDF_D0_RK5(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogCDF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: logCDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowetoLogCDF_D0_RK4(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogCDF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: logCDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowetoLogCDF_D0_RK3(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogCDF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: logCDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowetoLogCDF_D0_RK2(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogCDF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: logCDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowetoLogCDF_D0_RK1(logCDF, logx, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogCDF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: logCDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Inverse Cumulative Distribution (Quantile) Function of
    !>  the (Truncated) Poweto distribution for an input `logCDF` within the support of
    !>  the distribution \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[in]  logCDF      :   The input scalar (or array of the same shape as other array-like arguments)
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the CDF value
    !>                              within the support of the distribution at which the quantile must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `logCDF`,
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution,
    !>                              representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$. It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \return
    !>  `logx`                  :   The output scalar of the same type and kind the input argument `logCDF`,
    !>                              containing the natural logarithm of the quantile of the distribution at the specified CDF.<br>
    !>
    !>  \interface{getPowetoLogQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogQuan
    !>
    !>      logx = getPowetoLogQuan(logCDF, alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logCDF <= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPowetoLogQuan](@ref pm_distPoweto::setPowetoLogQuan)<br>
    !>
    !>  \example{getPowetoLogQuan}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogQuan/main.F90
    !>  \compilef{getPowetoLogQuan}
    !>  \output{getPowetoLogQuan}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogQuan/main.out.F90
    !>  \postproc{getPowetoLogQuan}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogQuan/main.py
    !>  \vis{getPowetoLogQuan}
    !>  \image html pm_distPoweto/getPowetoLogQuan/getPowetoLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPowetoLogQuan_D0_RK5(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogQuan_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logx
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPowetoLogQuan_D0_RK4(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogQuan_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logx
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPowetoLogQuan_D0_RK3(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogQuan_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logx
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPowetoLogQuan_D0_RK2(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogQuan_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logx
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPowetoLogQuan_D0_RK1(logCDF, alpha, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogQuan_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logMaxX
        real(RKC)                                           :: logx
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Inverse Cumulative Distribution (Quantile) Function of
    !>  the (Truncated) Poweto distribution for an input `logCDF` within the support of
    !>  the distribution \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[out] logx        :   The output scalar of the same type and kind the input argument `logCDF`, containing the natural logarithm of the CDF of the distribution.<br>
    !>  \param[in]  logCDF      :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  alpha       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of <br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same shape as other array-like arguments) of the same type, kind as `logx`,
    !>                              containing the natural logarithm of the normalization factor (\f$\eta\f$) of power-law function of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters significantly improves the runtime performance.<br>
    !>                              (**optional**, default = [getPowetoLogCDFNF(alpha, logMinX)](@ref pm_distPareto::getParetoLogCDFNF). It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \interface{getPowetoLogQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogQuan
    !>
    !>      call setPowetoLogQuan(logx, logCDF, alpha, logMinX = logMinX, logCDFNF = logCDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logCDFNF)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logCDF <= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowetoLogQuan](@ref pm_distPoweto::getPowetoLogQuan)<br>
    !>
    !>  \example{getPowetoLogQuan}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogQuan/main.F90
    !>  \compilef{getPowetoLogQuan}
    !>  \output{getPowetoLogQuan}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogQuan/main.out.F90
    !>  \postproc{getPowetoLogQuan}
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogQuan/main.py
    !>  \vis{getPowetoLogQuan}
    !>  \image html pm_distPoweto/getPowetoLogQuan/getPowetoLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowetoLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowetoLogQuan_D0_RK5(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogQuan_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: logx
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowetoLogQuan_D0_RK4(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogQuan_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: logx
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowetoLogQuan_D0_RK3(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogQuan_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: logx
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowetoLogQuan_D0_RK2(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogQuan_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: logx
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowetoLogQuan_D0_RK1(logx, logCDF, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogQuan_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: logx
        real(RKC)   , intent(in)                            :: logCDF
        real(RKC)   , intent(in)                            :: alpha
        real(RKC)   , intent(in)                , optional  :: logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of the natural logarithm(s) of random value(s) from the
    !>  <b>(Truncated) Poweto distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.<br>
    !>
    !>  \param[in]  alpha       :   The input scalar (or array of the same shape as other array-like arguments) of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the shape parameter of the distribution (i.e., the exponent of the power-law function).<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution,
    !>                              representing the minimum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$-\infty\f$. It must be present if \f$\alpha \leq 0\f$ holds.)
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution,
    !>                              representing the maximum of the support of the distribution.<br>
    !>                              (**optional**, default = \f$+\infty\f$. It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \return
    !>  `logRand`               :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `alpha`, containing the random value(s) from the specified distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogRand
    !>
    !>      logRand = getPowetoLogRand(alpha, logMinX = logMinX, logMaxX = logMaxX)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPowetoLogRand](@ref pm_distPoweto::setPowetoLogRand)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogRand/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distPoweto/getPowetoLogRand/main.py
    !>  \vis
    !>  \image html pm_distPoweto/getPowetoLogRand/getPowetoLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getPowetoLogRand_D0_RK5(alpha, logMinX, logMaxX, logCDFNF) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogRand_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX, logCDFNF
        real(RKC)                                   :: logRand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getPowetoLogRand_D0_RK4(alpha, logMinX, logMaxX, logCDFNF) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogRand_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX, logCDFNF
        real(RKC)                                   :: logRand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getPowetoLogRand_D0_RK3(alpha, logMinX, logMaxX, logCDFNF) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogRand_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX, logCDFNF
        real(RKC)                                   :: logRand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getPowetoLogRand_D0_RK2(alpha, logMinX, logMaxX, logCDFNF) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogRand_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX, logCDFNF
        real(RKC)                                   :: logRand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getPowetoLogRand_D0_RK1(alpha, logMinX, logMaxX, logCDFNF) result(logRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogRand_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logMaxX, logCDFNF
        real(RKC)                                   :: logRand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the natural logarithm(s) of random value(s) from the
    !>  <b>(Truncated) Poweto distribution</b> with parameters \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more
    !>  information on generating random numbers from the (Truncated) Poweto distribution.
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
    !>                              (**optional**, default = \f$-\infty\f$. It can be present <b>if and only if</b> `logCDFNF` is also present.)
    !>  \param[in]  logCDFNF    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `alpha`,
    !>                              containing the natural logarithm of the normalization factor of the CDF of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getPowetoLogCDFNF(alpha, logMinX, logMaxX)](@ref pm_distPoweto::getPowetoLogCDFNF).<br>
    !>                              (**optional**, default = [getPowetoLogCDFNF(alpha, logMinX)](@ref pm_distPareto::getParetoLogCDFNF). It must be present if \f$0 \leq \alpha\f$ holds.)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: setPowetoLogRand
    !>
    !>      call setPowetoLogRand(logRand, negExpRand, alpha, logMinX = logMinX, logCDFNF = logCDFNF)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < alpha .or. present(logMinX)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha < 0 .or. present(logMaxX)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  The conditions `negExpRand <= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPowetoLogRand](@ref pm_distPoweto::getPowetoLogRand)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogRand/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distPoweto/setPowetoLogRand/main.py
    !>  \vis
    !>  \image html pm_distPoweto/setPowetoLogRand/setPowetoLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow
    !>  This interface can be extended to support vector-like `logRand` arguments other than the `elemental` approach.<br>
    !>  Such an extension would be sensible only if the new interface improves the performance against the `elemental` approach.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowetoLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPowetoLogRand_D0_RK5(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogRand_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logCDFNF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPowetoLogRand_D0_RK4(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogRand_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logCDFNF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPowetoLogRand_D0_RK3(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogRand_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logCDFNF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPowetoLogRand_D0_RK2(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogRand_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logCDFNF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPowetoLogRand_D0_RK1(logRand, negExpRand, alpha, logMinX, logCDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogRand_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: negExpRand, alpha
        real(RKC)   , intent(in)    , optional      :: logMinX, logCDFNF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distPoweto
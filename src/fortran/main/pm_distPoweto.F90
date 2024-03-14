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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>(Truncated) Power/Pareto distribution</b> (hence the name **Poweto**).
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>(Truncated) Power/Pareto distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the **Cumulative Distribution Function (**CDF**)
!>      <li>    the random number generation from the distribution (**RNG**)
!>      <li>    the **Inverse Cumulative Distribution Function (ICDF)** or the **Quantile Function**
!>  </ol>
!>
!>  The (Truncated) Poweto distribution is better known by its special case: the <b>Broken Power Law distribution</b>.<br>
!>
!>  The **PDF** of an \f$n\f$-piece (Truncated) Poweto distribution (with \f$n > 0\f$) over a strictly-positive support
!>  \f$x \in [0 < x_\mathrm{lim,1}, x_\mathrm{lim,n+1} \leq +\infty)\f$ is defined with the **ascending-ordered** scale
!>  vector \f$x_\mathrm{lim} = \{ x_\mathrm{lim,i}, i = 1 : n + 1\}\f$ and shape vector \f$\alpha = \{\alpha_i, i = 1 : n \}\f$ as,
!>  \f{equation}{
!>      \large
!>      \pi(x | \alpha, x_\mathrm{lim}) =
!>      \begin{cases}
!>          \eta_1 x^{\alpha_1 - 1} &, ~ 0 \leq x_\mathrm{lim,1}  \leq x < x_\mathrm{lim,2} ~, \\
!>          \ldots \\
!>          \eta_i x^{\alpha_i - 1} &, ~        x_\mathrm{lim,i}  \leq x < x_\mathrm{lim,i+1} ~, \\
!>          \ldots \\
!>          \eta_n x^{\alpha_n - 1} &, ~        x_\mathrm{lim,n}  \leq x < x_\mathrm{lim,n+1} \leq +\infty ~,
!>      \end{cases}
!>  \f}
!>  where the conditions \f$\{\alpha_i \in \mathbb{R}, i = 1 : n\}\f$ and \f$\alpha_1 > 0 \lor x_\mathrm{lim,1} > 0\f$ and \f$\alpha_n < 0 \lor x_\mathrm{lim,n+1} < +\infty\f$ must hold.<br>
!>  These conditions must be met for the PDF to be normalizable.<br>
!>
!>  The component normalization factors \f$\eta_i\f$ are computed according to the following relationship,<br>
!>  \f{equation}{
!>      \large
!>      \eta_i = \eta_1 \prod_{j = 2}^i ~ x_\mathrm{lim,j}^{(\alpha_{j - 1} - \alpha_j)} ~, ~i = 2:n
!>  \f}
!>  where \f$\eta_1\f$ is a normalization factor that properly normalizes the integral of the PDF over its support to unity.<br>
!>
!>  The corresponding **Cumulative Distribution Function (CDF)** of the (Truncated) Poweto is,<br>
!>
!>  \f{equation}{
!>      \large
!>      \mathrm{CDF}(x | \alpha, x_\mathrm{lim}) = \sum_{i = 1}^{n : x_\mathrm{min,n} < x} ~ S_i ~,
!>  \f}
!>
!>  where \f$S_i ~:~ 1 \leq i \leq n\f$ is expressed via the following,<br>
!>
!>  \f{equation}{
!>      \large
!>      S_i =
!>      \begin{cases}
!>          \frac{\eta_i}{\alpha_i} \log\big( \frac{\min(x, x_{i+1})}{x_i} \big)                &, ~\alpha_i = 0 ~, \\
!>          \frac{\eta_i}{\alpha_i} \big( \min(x, x_{i+1})^{\alpha_i} - x_{i}^{\alpha_i} \big)  &, ~\alpha_i \neq 0 ~,
!>      \end{cases}
!>  \f}
!>
!>  where \f$\eta_i\f$ are the [normalization factors](@ref pm_distPoweto::getPowetoLogPDFNF) of the [Power-Law components](@ref pm_distPoweto) of the distribution.<br>
!>
!>  \see
!>  [pm_distPoweto](@ref pm_distPoweto)<br>
!>  [pm_distPoweto](@ref pm_distPoweto)<br>
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
    !>  Generate and return the natural logarithm of the normalization factors of the components of the Probability Density Function (PDF)
    !>  of the (Truncated) Poweto distribution for the input parameter vectors \f$(\alpha, x_\mathrm{lim})\f$.
    !>
    !>  \brief
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for the definition of the normalization factors.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factors of the (Truncated) Poweto distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the (Truncated) Poweto PDF to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  Alpha       :   The input vector of type `real` of kind \RKALL, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter(s) of the distribution (i.e., the exponents (plus one) of the power-law components of the distribution).<br>
    !>  \param[in]  LogLimX     :   The input vector of the same type and kind as `Alpha`, of size `size(Alpha) + 1` containing the natural logarithm of
    !>                              the scale parameters (i.e., the break points, or the limits) of the `n` power-law components of the distribution.<br>
    !>                              The scale parameter(s) must be in <b>ascending order</b>, such that `LogLimX(1) <= x <= LogLimX(size(LogLimX))`.<br>
    !>                              Setting `LogLimX(1) <= -log(huge(LogLimX))` effectively implies a left-opened semi-infinite support for the distribution.<br>
    !>                              Setting `LogLimX(size(LogLimX)) >= log(huge(LogLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>  \param[out] CumSumArea  :   The output vector of the same type, kind, and size as `LogLimX`, each element of which corresponds to the cumulative
    !>                              area underneath the distribution from the minimum of the support `exp(LogLimX(1))` to the corresponding element of `exp(LogLimX)`.<br>
    !>                              By definition, the conditions `CumSumArea(1) == 0.` and `CumSumArea(size(CumSumArea)) == 1.` and [isAscending(CumSumArea)](@ref pm_arraySort::isAscending) hold.<br>
    !>                              This output vector is a side-product of the computation of the normalization factors.<br>
    !>                              It is also required for [random number generation from the (Truncated) Poweto distribution](@ref DistPowetoRand_pmod).<br>
    !>                              Precomputing and supplying this vector to the random number generator routines significantly improves the runtime performance.<br>
    !>                              (**optional**. If missing, it will be computed implicitly within the algorithm and discarded upon return.)
    !>
    !>  \return
    !>  `LogNormFac`            :   The output vector of the same type, kind, and size as the input argument `Alpha`,
    !>                              containing the natural logarithm of the normalization factors of the power-law components of the (Truncated) Poweto distribution.<br>
    !>
    !>  \interface{getPowetoLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoLogPDFNF
    !>
    !>      LogNormFac(1:n) = getPowetoLogPDFNF(Alpha(1:n), LogLimX(1:n+1)) ! Poweto truncated at the upper limit `exp(LogLimX(n+1))`.
    !>      LogNormFac(1:n) = getPowetoLogPDFNF(Alpha(1:n), LogLimX(1:n+1), CumSumArea(1:n+1)) ! Poweto truncated at the upper limit `exp(LogLimX(n+1))`.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(Alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogLimX) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(CumSumArea) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogNormFac) == size(Alpha)` must hold for the corresponding input arguments.<br>
    !>  The condition `Alpha(1) > 0 .or. LogLimX(1) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `Alpha(size(Alpha)) < 0 .or. LogLimX(size(LogLimX)) < log(huge(LogLimX))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the output argument `CumSumArea` is present.<br>
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
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPowetoLogPDFNFALD_RK5(Alpha, LogLimX) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPowetoLogPDFNFALD_RK4(Alpha, LogLimX) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPowetoLogPDFNFALD_RK3(Alpha, LogLimX) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPowetoLogPDFNFALD_RK2(Alpha, LogLimX) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPowetoLogPDFNFALD_RK1(Alpha, LogLimX) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPowetoLogPDFNFALC_RK5(Alpha, LogLimX, CumSumArea) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALC_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)   , intent(out)   , contiguous    :: CumSumArea(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getPowetoLogPDFNFALC_RK4(Alpha, LogLimX, CumSumArea) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALC_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)   , intent(out)   , contiguous    :: CumSumArea(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getPowetoLogPDFNFALC_RK3(Alpha, LogLimX, CumSumArea) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALC_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)   , intent(out)   , contiguous    :: CumSumArea(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getPowetoLogPDFNFALC_RK2(Alpha, LogLimX, CumSumArea) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALC_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)   , intent(out)   , contiguous    :: CumSumArea(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getPowetoLogPDFNFALC_RK1(Alpha, LogLimX, CumSumArea) result(LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDFNFALC_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    , contiguous    :: Alpha(:), LogLimX(:)
        real(RKC)   , intent(out)   , contiguous    :: CumSumArea(:)
        real(RKC)                                   :: LogNormFac(size(Alpha, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the (Truncated) Poweto distribution for an input `logx`
    !>  within the support of the distribution `LogLimX(1) <= logx <= LogLimX(n+1)`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the `x` value within the support of the distribution at which the PDF must be computed.<br>
    !>  \param[in]  Alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  LogLimX     :   The input vector of the same type and kind as `logx`, of size `size(Alpha) + 1` containing the natural logarithm of the
    !>                              limits of the `n` power-law components of the distribution in ascending order, such that `LogLimX(1) <= x <= LogLimX(size(LogLimX))`.<br>
    !>                              Setting `LogLimX(size(LogLimX)) >= log(huge(LogLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>                              containing the natural logarithm of the scale parameters (i.e., the break points, or the minimum values of the power-law components) of the distribution.<br>
    !>  \param[in]  LogNormFac  :   The input vector of the same type, kind, and size as `Alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.<br>
    !>                              (**optional**, default = [getPowetoLogPDFNF(Alpha, LogLimX)](@ref pm_distPoweto::getPowetoLogPDFNF))
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
    !>      logPDF = getPowetoLogPDF(logx, Alpha(1:n), LogLimX(1:n+1), LogNormFac = LogNormFac(1:n))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(Alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogLimX) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogNormFac) == size(Alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `LogLimX(1) <= logx .and. logx < LogLimX(size(LogLimX))` must hold for the corresponding input arguments.<br>
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
    PURE module function getPowetoLogPDF_RK5(logx, Alpha, LogLimX, LogNormFac) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous, optional  :: LogNormFac(:)
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getPowetoLogPDF_RK4(logx, Alpha, LogLimX, LogNormFac) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous, optional  :: LogNormFac(:)
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getPowetoLogPDF_RK3(logx, Alpha, LogLimX, LogNormFac) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous, optional  :: LogNormFac(:)
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getPowetoLogPDF_RK2(logx, Alpha, LogLimX, LogNormFac) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous, optional  :: LogNormFac(:)
        real(RKC)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getPowetoLogPDF_RK1(logx, Alpha, LogLimX, LogNormFac) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoLogPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous, optional  :: LogNormFac(:)
        real(RKC)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the (Truncated) Poweto distribution for an input `logx`
    !>  within the support of the distribution `LogLimX(1) <= logx <= LogLimX(n+1)`.
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar of the same type and kind the input argument `logx`, containing the natural logarithm of the PDF of the distribution.<br>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the point at which the PDF must be computed.<br>
    !>  \param[in]  Alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  LogLimX     :   The input vector of the same type and kind as `logx`, of size `size(Alpha) + 1` containing the natural logarithm of the
    !>                              limits of the `n` power-law components of the distribution in ascending order, such that `LogLimX(1) <= x <= LogLimX(size(LogLimX))`.<br>
    !>                              Setting `LogLimX(size(LogLimX)) >= log(huge(LogLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>                              containing the natural logarithm of the scale parameters (i.e., the break points, or the minimum values of the power-law components) of the distribution.<br>
    !>  \param[in]  LogNormFac  :   The input vector of the same type, kind, and size as `Alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.
    !>
    !>  \interface{setPowetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: setPowetoLogPDF
    !>
    !>      call setPowetoLogPDF(logPDF, logx, Alpha(1:n), bin            , LogNormFac(1:n))
    !>      call setPowetoLogPDF(logPDF, logx, Alpha(1:n), LogLimX(1:n+1) , LogNormFac(1:n))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(Alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogLimX) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogNormFac) == size(Alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `LogLimX(1) <= logx .and. logx < LogLimX(size(LogLimX))` must hold for the corresponding input arguments.<br>
    !>  The conditions `1 <= bin .and. bin <= size(Alpha)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPowetoLogPDF](@ref pm_distPoweto::setPowetoLogPDF)<br>
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
    PURE module subroutine setPowetoLogPDFALL_D0_RK5(logPDF, logx, Alpha, LogLimX, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFALL_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPowetoLogPDFALL_D0_RK4(logPDF, logx, Alpha, LogLimX, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFALL_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPowetoLogPDFALL_D0_RK3(logPDF, logx, Alpha, LogLimX, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFALL_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPowetoLogPDFALL_D0_RK2(logPDF, logx, Alpha, LogLimX, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFALL_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPowetoLogPDFALL_D0_RK1(logPDF, logx, Alpha, LogLimX, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFALL_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPowetoLogPDFBAN_D0_RK5(logPDF, logx, Alpha, bin, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFBAN_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPowetoLogPDFBAN_D0_RK4(logPDF, logx, Alpha, bin, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFBAN_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPowetoLogPDFBAN_D0_RK3(logPDF, logx, Alpha, bin, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFBAN_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPowetoLogPDFBAN_D0_RK2(logPDF, logx, Alpha, bin, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFBAN_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPowetoLogPDFBAN_D0_RK1(logPDF, logx, Alpha, bin, LogNormFac)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoLogPDFBAN_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: logPDF
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the (Truncated) Poweto distribution for an input `logx`
    !>  within the support of the distribution `LogLimX(1) <= logx <= LogLimX(size(LogLimX))`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the `x` value within the support of the distribution at which the CDF must be computed.<br>
    !>  \param[in]  Alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  LogLimX     :   The input vector of the same type and kind as `Alpha`, of size `size(Alpha) + 1` containing the natural logarithm of
    !>                              the scale parameters (i.e., the break points, or the limits) of the `n` power-law components of the distribution.<br>
    !>                              The scale parameter(s) must be in <b>ascending order</b>, such that `LogLimX(1) <= x <= LogLimX(size(LogLimX))`.<br>
    !>                              Setting `LogLimX(1) <= -log(huge(LogLimX))` effectively implies a left-opened semi-infinite support for the distribution.<br>
    !>                              Setting `LogLimX(size(LogLimX)) >= log(huge(LogLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>  \param[in]  LogNormFac  :   The input vector of the same type, kind, and size as `Alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.<br>
    !>                              (**optional**. It must be present <b>if and only if</b> `CumSumArea` is also present. Default = [getPowetoLogPDFNF(Alpha, LogLimX)](@ref pm_distPoweto::getPowetoLogPDFNF))
    !>  \param[in]  CumSumArea  :   The output vector of the same type, kind, and size as `LogLimX`, each element of which corresponds to cumulative
    !>                              area underneath the distribution from the minimum of the support `exp(LogLimX(1))` to the corresponding element of `exp(LogLimX)`.<br>
    !>                              By definition, `CumSumArea(1) == 0.` and `CumSumArea(size(CumSumArea)) == 1.`, and [isAscending(CumSumArea)](@ref pm_arraySort::isAscending) hold.<br>
    !>                              This output vector is a side-product of the computation of the normalization factors.<br>
    !>                              It is also required for [random number generation from the (Truncated) Poweto distribution](@ref DistPowetoRand_pmod).<br>
    !>                              Precomputing and supplying this vector to the random number generator routines significantly improves the runtime performance.<br>
    !>                              (**optional**. It must be present <b>if and only if</b> `LogNormFac` is also present. 
    !>                              The default is set by [getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea)](@ref pm_distPoweto::getPowetoLogPDFNF).)<br>
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar of the same type and kind the input argument `logx`,
    !>                              containing the Cumulative Distribution Function (CDF) at the specified point within the support of the distribution.<br>
    !>
    !>  \interface{getPowetoCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: getPowetoCDF
    !>
    !>      cdf = getPowetoCDF(logx, Alpha(1:n), LogLimX(1:n+1))
    !>      cdf = getPowetoCDF(logx, Alpha(1:n), LogLimX(1:n+1), LogNormFac(1:n), CumSumArea(1:n+1))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(Alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogLimX) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(CumSumArea) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogNormFac) == size(Alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `LogLimX(1) <= logx .and. logx <= LogLimX(size(LogLimX))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setPowetoCDF](@ref pm_distPoweto::setPowetoCDF)<br>
    !>
    !>  \example{getPowetoCDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoCDF/main.F90
    !>  \compilef{getPowetoCDF}
    !>  \output{getPowetoCDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoCDF/main.out.F90
    !>  \postproc{getPowetoCDF}
    !>  \include{lineno} example/pm_distPoweto/getPowetoCDF/main.py
    !>  \vis{getPowetoCDF}
    !>  \image html pm_distPoweto/getPowetoCDF/getPowetoCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{getPowetoCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPowetoCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getPowetoCDFALDD_RK5(logx, Alpha, LogLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALDD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK4_ENABLED
    impure module function getPowetoCDFALDD_RK4(logx, Alpha, LogLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALDD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK3_ENABLED
    impure module function getPowetoCDFALDD_RK3(logx, Alpha, LogLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALDD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK2_ENABLED
    impure module function getPowetoCDFALDD_RK2(logx, Alpha, LogLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALDD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK1_ENABLED
    impure module function getPowetoCDFALDD_RK1(logx, Alpha, LogLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALDD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)                                           :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getPowetoCDFALLC_RK5(logx, Alpha, LogLimX, LogNormFac, CumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALLC_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK4_ENABLED
    impure module function getPowetoCDFALLC_RK4(logx, Alpha, LogLimX, LogNormFac, CumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALLC_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK3_ENABLED
    impure module function getPowetoCDFALLC_RK3(logx, Alpha, LogLimX, LogNormFac, CumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALLC_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK2_ENABLED
    impure module function getPowetoCDFALLC_RK2(logx, Alpha, LogLimX, LogNormFac, CumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALLC_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        real(RKC)                                           :: cdf
    end function
#endif

#if RK1_ENABLED
    impure module function getPowetoCDFALLC_RK1(logx, Alpha, LogLimX, LogNormFac, CumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPowetoCDFALLC_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        real(RKC)                                           :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the (Truncated) Poweto distribution for an input `logx`
    !>  within the support of the distribution `LogLimX(1) <= logx <= LogLimX(size(LogLimX))`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPoweto](@ref pm_distPoweto) for more information on the (Truncated) Poweto distribution.<br>
    !>
    !>  \param[out] cdf         :   The output scalar of the same type and kind the input argument `logx`, containing the Cumulative Distribution Function (CDF) at the specified `logx`.<br>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the point at which the CDF must be computed.<br>
    !>  \param[in]  Alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  LogLimX     :   The input vector of the same type and kind as `Alpha`, of size `size(Alpha) + 1` containing the natural logarithm of
    !>                              the scale parameters (i.e., the break points, or the limits) of the `n` power-law components of the distribution.<br>
    !>                              The scale parameter(s) must be in <b>ascending order</b>, such that `LogLimX(1) <= x <= LogLimX(size(LogLimX))`.<br>
    !>                              Setting `LogLimX(1) <= -log(huge(LogLimX))` effectively implies a left-opened semi-infinite support for the distribution.<br>
    !>                              Setting `LogLimX(size(LogLimX)) >= log(huge(LogLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>  \param[in]  LogNormFac  :   The input vector of the same type, kind, and size as `Alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) Poweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.<br>
    !>                              This input vector can be readily obtained by calling [getPowetoLogPDFNF(Alpha, LogLimX)](@ref pm_distPoweto::getPowetoLogPDFNF).<br>
    !>  \param[in]  CumSumArea  :   The output vector of the same type, kind, and size as `LogLimX`, each element of which corresponds to cumulative
    !>                              area underneath the distribution from the minimum of the support `exp(LogLimX(1))` to the corresponding element of `exp(LogLimX)`.<br>
    !>                              By definition, `CumSumArea(1) == 0.` and `CumSumArea(size(CumSumArea)) == 1.`, and [isAscending(CumSumArea)](@ref pm_arraySort::isAscending) hold.<br>
    !>                              This output vector is a side-product of the computation of the normalization factors.<br>
    !>                              It is also required for [random number generation from the (Truncated) Poweto distribution](@ref DistPowetoRand_pmod).<br>
    !>                              Precomputing and supplying this vector to the random number generator routines significantly improves the runtime performance.<br>
    !>                              This vector can be readily obtained by calling [getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea)](@ref pm_distPoweto::getPowetoLogPDFNF).<br>
    !>
    !>  \interface{setPowetoCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPoweto, only: setPowetoCDF
    !>
    !>      call setPowetoCDF(cdf, logx, Alpha(1:n), LogLimX(1:n+1), LogNormFac(1:n), CumSumArea(1:n+1))
    !>      call setPowetoCDF(cdf, logx, Alpha(1:n), LogLimX(1:n+1), LogNormFac(1:n), CumSumArea(1:n+1), bin)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(Alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogLimX) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(CumSumArea) == size(Alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(LogNormFac) == size(Alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `LogLimX(1) <= logx .and. logx <= LogLimX(size(LogLimX))` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < bin .and. bin < size(LogLimX)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPowetoCDF](@ref pm_distPoweto::getPowetoCDF)<br>
    !>
    !>  \example{setPowetoCDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoCDF/main.F90
    !>  \compilef{setPowetoCDF}
    !>  \output{setPowetoCDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoCDF/main.out.F90
    !>  \postproc{setPowetoCDF}
    !>  \include{lineno} example/pm_distPoweto/setPowetoCDF/main.py
    !>  \vis{setPowetoCDF}
    !>  \image html pm_distPoweto/setPowetoCDF/setPowetoCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPoweto](@ref test_pm_distPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setPowetoCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPowetoCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPowetoCDFMAN_D0_RK5(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFMAN_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPowetoCDFMAN_D0_RK4(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFMAN_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPowetoCDFMAN_D0_RK3(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFMAN_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPowetoCDFMAN_D0_RK2(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFMAN_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPowetoCDFMAN_D0_RK1(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFMAN_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPowetoCDFBAN_D0_RK5(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFBAN_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPowetoCDFBAN_D0_RK4(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFBAN_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPowetoCDFBAN_D0_RK3(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFBAN_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPowetoCDFBAN_D0_RK2(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFBAN_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPowetoCDFBAN_D0_RK1(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPowetoCDFBAN_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                           :: cdf
        real(RKC)   , intent(in)                            :: logx
        real(RKC)   , intent(in)    , contiguous            :: Alpha(:)
        real(RKC)   , intent(in)    , contiguous            :: LogLimX(:)
        real(RKC)   , intent(in)    , contiguous            :: LogNormFac(:)
        real(RKC)   , intent(in)    , contiguous            :: CumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distPoweto
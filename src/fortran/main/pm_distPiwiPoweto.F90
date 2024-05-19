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
!>  related to the <b>(Truncated) PieceWise Power/Pareto distribution</b> (hence the name **PiwiPoweto**).
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>(Truncated) PieceWise Power/Pareto distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The (Truncated) PiwiPoweto distribution is better known by its special case: the <b>Broken Power Law distribution</b>.<br>
!>
!>  The **PDF** of an \f$n\f$-piece (Truncated) PiwiPoweto distribution (with \f$n > 0\f$) over a strictly-positive support
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
!>  where the conditions \f$\{\alpha_i \in \mathbb{R}, i = 1 : n\}\f$ and \f$\alpha_1 > 0 \lor x_\mathrm{lim,1} > 0\f$
!>  and \f$\alpha_n < 0 \lor x_\mathrm{lim,n+1} < +\infty\f$ must hold.<br>
!>  These conditions must be met for the PDF to be normalizable.<br>
!>
!>  The component normalization factors \f$\eta_i\f$ are computed according to the following relationship,<br>
!>  \f{equation}{
!>      \large
!>      \eta_i = \eta_1 \prod_{j = 2}^i ~ x_\mathrm{lim,j}^{(\alpha_{j - 1} - \alpha_j)} ~, ~i = 2:n
!>  \f}
!>  where \f$\eta_1\f$ is a normalization factor that properly normalizes the integral of the PDF over its support to unity.<br>
!>
!>  The corresponding **Cumulative Distribution Function (CDF)** of the (Truncated) PiwiPoweto is,<br>
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
!>  where \f$\eta_i\f$ are the [normalization factors](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDFNF) of the [Power-Law components](@ref pm_distPiwiPoweto) of the distribution.<br>
!>
!>  \see
!>  [pm_distPower](@ref pm_distPower)<br>
!>  [pm_distPareto](@ref pm_distPareto)<br>
!>
!>  \test
!>  [test_pm_distPiwiPoweto](@ref test_pm_distPiwiPoweto)
!>
!>  \todo
!>  Generic interfaces for computing the logarithm of CDF robustly (without numerical rounding) must be added in the future.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distPiwiPoweto

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distPiwiPoweto"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization factors of the components of the Probability Density Function (PDF)
    !>  of the (Truncated) PiwiPoweto distribution for the input parameter vectors \f$(\alpha, x_\mathrm{lim})\f$.
    !>
    !>  \brief
    !>  See the documentation of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto) for the definition of the normalization factors.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factors of the (Truncated) PiwiPoweto distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the (Truncated) PiwiPoweto PDF to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  alpha       :   The input vector of type `real` of kind \RKALL, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter(s) of the distribution (i.e., the exponents (plus one) of the power-law components of the distribution).<br>
    !>  \param[in]  logLimX     :   The input vector of the same type and kind as `alpha`, of size `size(alpha) + 1` containing the natural logarithm of
    !>                              the scale parameters (i.e., the break points, or the limits) of the `n` power-law components of the distribution.<br>
    !>                              The scale parameter(s) must be in <b>ascending order</b>, such that `logLimX(1) <= x <= logLimX(size(logLimX))`.<br>
    !>                              Setting `logLimX(1) <= -log(huge(logLimX))` effectively implies a left-opened semi-infinite support for the distribution.<br>
    !>                              Setting `logLimX(size(logLimX)) >= log(huge(logLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>  \param[out] cumSumArea  :   The output vector of the same type, kind, and size as `logLimX`, each element of which corresponds to the cumulative
    !>                              area underneath the distribution from the minimum of the support `exp(logLimX(1))` to the corresponding element of `exp(logLimX)`.<br>
    !>                              By definition, the conditions `cumSumArea(1) == 0.` and `cumSumArea(size(cumSumArea)) == 1.` and [isAscending(cumSumArea)](@ref pm_arraySort::isAscending) hold.<br>
    !>                              This output vector is a side-product of the computation of the normalization factors.<br>
    !>                              It is also required for random number generation from the (Truncated) PiwiPoweto distribution.<br>
    !>                              Precomputing and supplying this vector to the random number generator routines significantly improves the runtime performance.<br>
    !>                              (**optional**. If missing, it will be computed implicitly within the algorithm and discarded upon return.)
    !>
    !>  \return
    !>  `logPDFNF`              :   The output vector of the same type, kind, and size as the input argument `alpha`,
    !>                              containing the natural logarithm of the normalization factors of the power-law components of the (Truncated) PiwiPoweto distribution.<br>
    !>
    !>  \interface{getPiwiPowetoLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distPiwiPoweto, only: getPiwiPowetoLogPDFNF
    !>
    !>      logPDFNF(1:n) = getPiwiPowetoLogPDFNF(alpha(1:n), logLimX(1:n+1)) ! PiwiPoweto truncated at the upper limit `exp(logLimX(n+1))`.
    !>      logPDFNF(1:n) = getPiwiPowetoLogPDFNF(alpha(1:n), logLimX(1:n+1), cumSumArea(1:n+1)) ! PiwiPoweto truncated at the upper limit `exp(logLimX(n+1))`.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logLimX) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(cumSumArea) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logPDFNF) == size(alpha)` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha(1) > 0 .or. logLimX(1) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha(size(alpha)) < 0 .or. logLimX(size(logLimX)) < log(huge(logLimX))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the output argument `cumSumArea` is present.<br>
    !>
    !>  \see
    !>  [getPiwiPowetoLogPDF](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDF)<br>
    !>  [setPiwiPowetoLogPDF](@ref pm_distPiwiPoweto::setPiwiPowetoLogPDF)<br>
    !>
    !>  \example{getPiwiPowetoLogPDFNF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoLogPDFNF/main.F90
    !>  \compilef{getPiwiPowetoLogPDFNF}
    !>  \output{getPiwiPowetoLogPDFNF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoLogPDFNF/main.out.F90
    !>  \postproc{getPiwiPowetoLogPDFNF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoLogPDFNF/main.py
    !>  \vis
    !>  \image html pm_distPiwiPoweto/getPiwiPowetoLogPDFNF/getPiwiPowetoLogPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPiwiPoweto](@ref test_pm_distPiwiPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getPiwiPowetoLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPiwiPowetoLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPiwiPowetoLogPDFNFALD_RK5(alpha, logLimX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPiwiPowetoLogPDFNFALD_RK4(alpha, logLimX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPiwiPowetoLogPDFNFALD_RK3(alpha, logLimX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPiwiPowetoLogPDFNFALD_RK2(alpha, logLimX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPiwiPowetoLogPDFNFALD_RK1(alpha, logLimX) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPiwiPowetoLogPDFNFALC_RK5(alpha, logLimX, cumSumArea) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALC_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)   , intent(out)   , contiguous    :: cumSumArea(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getPiwiPowetoLogPDFNFALC_RK4(alpha, logLimX, cumSumArea) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALC_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)   , intent(out)   , contiguous    :: cumSumArea(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getPiwiPowetoLogPDFNFALC_RK3(alpha, logLimX, cumSumArea) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALC_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)   , intent(out)   , contiguous    :: cumSumArea(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getPiwiPowetoLogPDFNFALC_RK2(alpha, logLimX, cumSumArea) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALC_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)   , intent(out)   , contiguous    :: cumSumArea(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getPiwiPowetoLogPDFNFALC_RK1(alpha, logLimX, cumSumArea) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDFNFALC_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: alpha(:), logLimX(:)
        real(RKG)   , intent(out)   , contiguous    :: cumSumArea(:)
        real(RKG)                                   :: logPDFNF(size(alpha, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the (Truncated) PiwiPoweto distribution for an input `logx`
    !>  within the support of the distribution `logLimX(1) <= logx <= logLimX(n+1)`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto) for more information on the (Truncated) PiwiPoweto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the `x` value within the support of the distribution at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  logLimX     :   The input vector of the same type and kind as `logx`, of size `size(alpha) + 1` containing the natural logarithm of the
    !>                              limits of the `n` power-law components of the distribution in ascending order, such that `logLimX(1) <= x <= logLimX(size(logLimX))`.<br>
    !>                              Setting `logLimX(size(logLimX)) >= log(huge(logLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>                              containing the natural logarithm of the scale parameters (i.e., the break points, or the minimum values of the power-law components) of the distribution.<br>
    !>  \param[in]  logPDFNF    :   The input vector of the same type, kind, and size as `alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) PiwiPoweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.<br>
    !>                              (**optional**, default = [getPiwiPowetoLogPDFNF(alpha, logLimX)](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDFNF))
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar of the same type and kind the input argument `logx`,
    !>                              containing the natural logarithm of the PDF of the distribution at the specified point within the support of the PDF.<br>
    !>
    !>  \interface{getPiwiPowetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPiwiPoweto, only: getPiwiPowetoLogPDF
    !>
    !>      logPDF = getPiwiPowetoLogPDF(logx, alpha(1:n), logLimX(1:n+1), logPDFNF = logPDFNF(1:n))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logLimX) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logPDFNF) == size(alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logLimX(1) <= logx .and. logx < logLimX(size(logLimX))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPiwiPowetoLogPDF](@ref pm_distPiwiPoweto::setPiwiPowetoLogPDF)<br>
    !>
    !>  \example{getPiwiPowetoLogPDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoLogPDF/main.F90
    !>  \compilef{getPiwiPowetoLogPDF}
    !>  \output{getPiwiPowetoLogPDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoLogPDF/main.out.F90
    !>  \postproc{getPiwiPowetoLogPDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoLogPDF/main.py
    !>  \vis{getPiwiPowetoLogPDF}
    !>  \image html pm_distPiwiPoweto/getPiwiPowetoLogPDF/getPiwiPowetoLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPiwiPoweto](@ref test_pm_distPiwiPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getPiwiPowetoLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPiwiPowetoLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPiwiPowetoLogPDF_RK5(logx, alpha, logLimX, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous, optional  :: logPDFNF(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getPiwiPowetoLogPDF_RK4(logx, alpha, logLimX, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous, optional  :: logPDFNF(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getPiwiPowetoLogPDF_RK3(logx, alpha, logLimX, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous, optional  :: logPDFNF(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getPiwiPowetoLogPDF_RK2(logx, alpha, logLimX, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous, optional  :: logPDFNF(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getPiwiPowetoLogPDF_RK1(logx, alpha, logLimX, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous, optional  :: logPDFNF(:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the (Truncated) PiwiPoweto distribution for an input `logx`
    !>  within the support of the distribution `logLimX(1) <= logx <= logLimX(n+1)`.
    !>
    !>  \details
    !>  See the documentation of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto) for more information on the (Truncated) PiwiPoweto distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar of the same type and kind the input argument `logx`, containing the natural logarithm of the PDF of the distribution.<br>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the point at which the PDF must be computed.<br>
    !>  \param[in]  alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  logLimX     :   The input vector of the same type and kind as `logx`, of size `size(alpha) + 1` containing the natural logarithm of the
    !>                              limits of the `n` power-law components of the distribution in ascending order, such that `logLimX(1) <= x <= logLimX(size(logLimX))`.<br>
    !>                              Setting `logLimX(size(logLimX)) >= log(huge(logLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>                              containing the natural logarithm of the scale parameters (i.e., the break points, or the minimum values of the power-law components) of the distribution.<br>
    !>  \param[in]  logPDFNF    :   The input vector of the same type, kind, and size as `alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) PiwiPoweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.
    !>
    !>  \interface{setPiwiPowetoLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPiwiPoweto, only: setPiwiPowetoLogPDF
    !>
    !>      call setPiwiPowetoLogPDF(logPDF, logx, alpha(1:n), bin            , logPDFNF(1:n))
    !>      call setPiwiPowetoLogPDF(logPDF, logx, alpha(1:n), logLimX(1:n+1) , logPDFNF(1:n))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logLimX) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logPDFNF) == size(alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logLimX(1) <= logx .and. logx < logLimX(size(logLimX))` must hold for the corresponding input arguments.<br>
    !>  The conditions `1 <= bin .and. bin <= size(alpha)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPiwiPowetoLogPDF](@ref pm_distPiwiPoweto::setPiwiPowetoLogPDF)<br>
    !>
    !>  \example{setPiwiPowetoLogPDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/setPiwiPowetoLogPDF/main.F90
    !>  \compilef{setPiwiPowetoLogPDF}
    !>  \output{setPiwiPowetoLogPDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/setPiwiPowetoLogPDF/main.out.F90
    !>  \postproc{setPiwiPowetoLogPDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/setPiwiPowetoLogPDF/main.py
    !>  \vis{setPiwiPowetoLogPDF}
    !>  \image html pm_distPiwiPoweto/setPiwiPowetoLogPDF/setPiwiPowetoLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPiwiPoweto](@ref test_pm_distPiwiPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setPiwiPowetoLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPiwiPowetoLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFALL_D0_RK5(logPDF, logx, alpha, logLimX, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFALL_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFALL_D0_RK4(logPDF, logx, alpha, logLimX, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFALL_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFALL_D0_RK3(logPDF, logx, alpha, logLimX, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFALL_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFALL_D0_RK2(logPDF, logx, alpha, logLimX, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFALL_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFALL_D0_RK1(logPDF, logx, alpha, logLimX, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFALL_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFBAN_D0_RK5(logPDF, logx, alpha, bin, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFBAN_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFBAN_D0_RK4(logPDF, logx, alpha, bin, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFBAN_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFBAN_D0_RK3(logPDF, logx, alpha, bin, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFBAN_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFBAN_D0_RK2(logPDF, logx, alpha, bin, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFBAN_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPiwiPowetoLogPDFBAN_D0_RK1(logPDF, logx, alpha, bin, logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoLogPDFBAN_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                           :: logPDF
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the (Truncated) PiwiPoweto distribution for an input `logx`
    !>  within the support of the distribution `logLimX(1) <= logx <= logLimX(size(logLimX))`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto) for more information on the (Truncated) PiwiPoweto distribution.
    !>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the `x` value within the support of the distribution at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  logLimX     :   The input vector of the same type and kind as `alpha`, of size `size(alpha) + 1` containing the natural logarithm of
    !>                              the scale parameters (i.e., the break points, or the limits) of the `n` power-law components of the distribution.<br>
    !>                              The scale parameter(s) must be in <b>ascending order</b>, such that `logLimX(1) <= x <= logLimX(size(logLimX))`.<br>
    !>                              Setting `logLimX(1) <= -log(huge(logLimX))` effectively implies a left-opened semi-infinite support for the distribution.<br>
    !>                              Setting `logLimX(size(logLimX)) >= log(huge(logLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>  \param[in]  logPDFNF    :   The input vector of the same type, kind, and size as `alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) PiwiPoweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.<br>
    !>                              (**optional**. It must be present <b>if and only if</b> `cumSumArea` is also present. Default = [getPiwiPowetoLogPDFNF(alpha, logLimX)](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDFNF))
    !>  \param[in]  cumSumArea  :   The output vector of the same type, kind, and size as `logLimX`, each element of which corresponds to cumulative
    !>                              area underneath the distribution from the minimum of the support `exp(logLimX(1))` to the corresponding element of `exp(logLimX)`.<br>
    !>                              By definition, `cumSumArea(1) == 0.` and `cumSumArea(size(cumSumArea)) == 1.`, and [isAscending(cumSumArea)](@ref pm_arraySort::isAscending) hold.<br>
    !>                              This output vector is a side-product of the computation of the normalization factors.<br>
    !>                              It is also required for random number generation from the (Truncated) PiwiPoweto distribution.<br>
    !>                              Precomputing and supplying this vector to the random number generator routines significantly improves the runtime performance.<br>
    !>                              (**optional**. It must be present <b>if and only if</b> `logPDFNF` is also present.
    !>                              The default is set by [getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea)](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDFNF).)<br>
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar of the same type and kind the input argument `logx`,
    !>                              containing the Cumulative Distribution Function (CDF) at the specified point within the support of the distribution.<br>
    !>
    !>  \interface{getPiwiPowetoCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPiwiPoweto, only: getPiwiPowetoCDF
    !>
    !>      cdf = getPiwiPowetoCDF(logx, alpha(1:n), logLimX(1:n+1))
    !>      cdf = getPiwiPowetoCDF(logx, alpha(1:n), logLimX(1:n+1), logPDFNF(1:n), cumSumArea(1:n+1))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logLimX) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(cumSumArea) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logPDFNF) == size(alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logLimX(1) <= logx .and. logx <= logLimX(size(logLimX))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setPiwiPowetoCDF](@ref pm_distPiwiPoweto::setPiwiPowetoCDF)<br>
    !>
    !>  \example{getPiwiPowetoCDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoCDF/main.F90
    !>  \compilef{getPiwiPowetoCDF}
    !>  \output{getPiwiPowetoCDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoCDF/main.out.F90
    !>  \postproc{getPiwiPowetoCDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/getPiwiPowetoCDF/main.py
    !>  \vis{getPiwiPowetoCDF}
    !>  \image html pm_distPiwiPoweto/getPiwiPowetoCDF/getPiwiPowetoCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPiwiPoweto](@ref test_pm_distPiwiPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getPiwiPowetoCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPiwiPowetoCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getPiwiPowetoCDFALDD_RK5(logx, alpha, logLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK4_ENABLED
    impure module function getPiwiPowetoCDFALDD_RK4(logx, alpha, logLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK3_ENABLED
    impure module function getPiwiPowetoCDFALDD_RK3(logx, alpha, logLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK2_ENABLED
    impure module function getPiwiPowetoCDFALDD_RK2(logx, alpha, logLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK1_ENABLED
    impure module function getPiwiPowetoCDFALDD_RK1(logx, alpha, logLimX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)                                           :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getPiwiPowetoCDFALLC_RK5(logx, alpha, logLimX, logPDFNF, cumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALLC_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK4_ENABLED
    impure module function getPiwiPowetoCDFALLC_RK4(logx, alpha, logLimX, logPDFNF, cumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALLC_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK3_ENABLED
    impure module function getPiwiPowetoCDFALLC_RK3(logx, alpha, logLimX, logPDFNF, cumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALLC_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK2_ENABLED
    impure module function getPiwiPowetoCDFALLC_RK2(logx, alpha, logLimX, logPDFNF, cumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALLC_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        real(RKG)                                           :: cdf
    end function
#endif

#if RK1_ENABLED
    impure module function getPiwiPowetoCDFALLC_RK1(logx, alpha, logLimX, logPDFNF, cumSumArea) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPiwiPowetoCDFALLC_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        real(RKG)                                           :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the (Truncated) PiwiPoweto distribution for an input `logx`
    !>  within the support of the distribution `logLimX(1) <= logx <= logLimX(size(logLimX))`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto) for more information on the (Truncated) PiwiPoweto distribution.<br>
    !>
    !>  \param[out] cdf         :   The output scalar of the same type and kind the input argument `logx`, containing the Cumulative Distribution Function (CDF) at the specified `logx`.<br>
    !>  \param[in]  logx        :   The input scalar of type `real` of kind \RKALL, containing the natural logarithm of the point at which the CDF must be computed.<br>
    !>  \param[in]  alpha       :   The input vector of the same type and kind as `logx`, of the same size `n` as the number of the power-law components of the distribution,
    !>                              containing the shape parameter of the distribution (i.e., the exponents of the power-law components of the distribution).<br>
    !>  \param[in]  logLimX     :   The input vector of the same type and kind as `alpha`, of size `size(alpha) + 1` containing the natural logarithm of
    !>                              the scale parameters (i.e., the break points, or the limits) of the `n` power-law components of the distribution.<br>
    !>                              The scale parameter(s) must be in <b>ascending order</b>, such that `logLimX(1) <= x <= logLimX(size(logLimX))`.<br>
    !>                              Setting `logLimX(1) <= -log(huge(logLimX))` effectively implies a left-opened semi-infinite support for the distribution.<br>
    !>                              Setting `logLimX(size(logLimX)) >= log(huge(logLimX))` effectively implies a right-opened semi-infinite support for the distribution.<br>
    !>  \param[in]  logPDFNF    :   The input vector of the same type, kind, and size as `alpha`, containing the natural logarithm of the normalization factors (\f$\eta\f$)
    !>                              of power-law components of the distribution of the (Truncated) PiwiPoweto distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(\alpha, x_\mathrm{lim})\f$ parameters
    !>                              significantly improves the runtime performance.<br>
    !>                              This input vector can be readily obtained by calling [getPiwiPowetoLogPDFNF(alpha, logLimX)](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDFNF).<br>
    !>  \param[in]  cumSumArea  :   The output vector of the same type, kind, and size as `logLimX`, each element of which corresponds to cumulative
    !>                              area underneath the distribution from the minimum of the support `exp(logLimX(1))` to the corresponding element of `exp(logLimX)`.<br>
    !>                              By definition, `cumSumArea(1) == 0.` and `cumSumArea(size(cumSumArea)) == 1.`, and [isAscending(cumSumArea)](@ref pm_arraySort::isAscending) hold.<br>
    !>                              This output vector is a side-product of the computation of the normalization factors.<br>
    !>                              It is also required for random number generation from the (Truncated) PiwiPoweto distribution.<br>
    !>                              Precomputing and supplying this vector to the random number generator routines significantly improves the runtime performance.<br>
    !>                              This vector can be readily obtained by calling [getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea)](@ref pm_distPiwiPoweto::getPiwiPowetoLogPDFNF).<br>
    !>
    !>  \interface{setPiwiPowetoCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPiwiPoweto, only: setPiwiPowetoCDF
    !>
    !>      call setPiwiPowetoCDF(cdf, logx, alpha(1:n), logLimX(1:n+1), logPDFNF(1:n), cumSumArea(1:n+1))
    !>      call setPiwiPowetoCDF(cdf, logx, alpha(1:n), logLimX(1:n+1), logPDFNF(1:n), cumSumArea(1:n+1), bin)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(alpha) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logLimX) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(cumSumArea) == size(alpha) + 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(logPDFNF) == size(alpha)` must hold for the corresponding input arguments.<br>
    !>  The conditions `logLimX(1) <= logx .and. logx <= logLimX(size(logLimX))` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < bin .and. bin < size(logLimX)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPiwiPowetoCDF](@ref pm_distPiwiPoweto::getPiwiPowetoCDF)<br>
    !>
    !>  \example{setPiwiPowetoCDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/setPiwiPowetoCDF/main.F90
    !>  \compilef{setPiwiPowetoCDF}
    !>  \output{setPiwiPowetoCDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/setPiwiPowetoCDF/main.out.F90
    !>  \postproc{setPiwiPowetoCDF}
    !>  \include{lineno} example/pm_distPiwiPoweto/setPiwiPowetoCDF/main.py
    !>  \vis{setPiwiPowetoCDF}
    !>  \image html pm_distPiwiPoweto/setPiwiPowetoCDF/setPiwiPowetoCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPiwiPoweto](@ref test_pm_distPiwiPoweto)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setPiwiPowetoCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPiwiPowetoCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPiwiPowetoCDFMAN_D0_RK5(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFMAN_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPiwiPowetoCDFMAN_D0_RK4(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFMAN_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPiwiPowetoCDFMAN_D0_RK3(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFMAN_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPiwiPowetoCDFMAN_D0_RK2(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFMAN_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPiwiPowetoCDFMAN_D0_RK1(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFMAN_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPiwiPowetoCDFBAN_D0_RK5(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFBAN_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPiwiPowetoCDFBAN_D0_RK4(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFBAN_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPiwiPowetoCDFBAN_D0_RK3(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFBAN_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPiwiPowetoCDFBAN_D0_RK2(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFBAN_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPiwiPowetoCDFBAN_D0_RK1(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea, bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPiwiPowetoCDFBAN_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                           :: cdf
        real(RKG)   , intent(in)                            :: logx
        real(RKG)   , intent(in)    , contiguous            :: alpha(:)
        real(RKG)   , intent(in)    , contiguous            :: logLimX(:)
        real(RKG)   , intent(in)    , contiguous            :: logPDFNF(:)
        real(RKG)   , intent(in)    , contiguous            :: cumSumArea(:)
        integer(IK) , intent(in)                            :: bin
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distPiwiPoweto
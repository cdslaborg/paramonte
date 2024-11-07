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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>univariate Normal distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>univariate Normal distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the Normal distribution is defined with the two location and scale parameters \f$(\mu \in (-\infty, +\infty), \sigma > 0)\f$ as,
!>  \f{equation}{
!>      \large
!>      \pi(x | \mu, \sigma) =
!>      \frac{1}{\sigma\sqrt{2\pi}}\exp\bigg( -\frac{\big(x - \mu\big)^2}{2\sigma^2} \bigg) ~,~ x \in (-\infty, +\infty) ~.
!>  \f}
!>
!>  The **CDF** of the Normal distribution is defined with the two location and scale parameters \f$(\mu \in (-\infty, +\infty), \sigma > 0)\f$ as,
!>  \f{equation}{
!>      \large
!>      \mathrm{CDF}(x | \mu, \sigma) =
!>      \frac{1}{2} \bigg[ 1 + \mathrm{erf} \bigg( \frac{x - \mu}{\sigma\sqrt{2}} \bigg) \bigg] ~,~ x \in (-\infty, +\infty) ~.
!>  \f}
!>
!>  **Quantile Function**<br>
!>
!>  In probability and statistics, the quantile function outputs the value of a
!>  random variable such that its probability is less than or equal to an input probability value.<br>
!>  Intuitively, the quantile function associates with a range at and below a probability input the
!>  likelihood that a random variable is realized in that range for some probability distribution.<br>
!>  It is also called the **percentile function**, **percent-point function** or **inverse cumulative distribution function (ICDF)**.<br>
!>
!>  See the documentation of [pm_distNorm](@ref pm_distNorm) for information on the CDF of the Normal distribution.<br>
!>
!>  The quantile function of the standard normal distribution is called the [**probit function**](https://en.wikipedia.org/wiki/Probit),
!>  and can be expressed in terms of the inverse error function:
!>  \f{equation}{
!>      \Phi^{-1}(p) = {\sqrt{2}}\ms{erf}^{-1}(2p - 1), \quad p \in (0,1) ~.
!>  \f}
!>  For a normal random variable with mean \f$\mu\f$ and variance \f$\sigma^{2}\f$, the quantile function is,
!>  \f{equation}{
!>      F^{-1}(p) = \mu + \sigma \Phi^{-1}(p) = \mu + \sigma \sqrt{2} \ms{erf}^{-1}(2p - 1), \quad p \in (0,1) ~.
!>  \f}
!>  The quantile \f$\Phi^{-1}(p)\f$ of the standard normal distribution is commonly denoted as \f$z_{p}\f$.<br>
!>  These values are used in hypothesis testing, construction of confidence intervals and **Q–Q** plots.<br>
!>  A normal random variable \f$X\f$ will exceed \f$\mu + z_{p}\sigma\f$ with probability \f$1 − p\f$,
!>  and will lie outside the interval \f$\mu \pm z_{p} \sigma\f$ with probability \f$2(1 - p)\f$.<br>
!>  In particular, the quantile \f$z_{0.975}\f$ is \f$1.96\f$.<br>
!>  Therefore a normal random variable will lie outside the interval \f$\mu \pm 1.96\sigma\f$ in only \f$5\%\f$ of cases.<br>
!>
!>  The following table gives the quantile \f$z_{p}\f$ such that \f$X\f$ will lie in the range \f$\mu \pm z_{p}\sigma\f$ with a specified probability \f$p\f$.<br>
!>  These values are useful to determine tolerance interval for sample averages and other statistical estimators with normal (or asymptotically normal) distributions.<br>
!>  Note that the following table shows \f$\sqrt{2} \ms{erf}^{-1}(p) = \Phi^{-1}\left(\frac{p + 1}{2}\right)\f$, not \f$\Phi^{-1}(p)\f$ as defined above.<br>
!>
!>  \f$p\f$     |   \f$z_{p}\f$
!>  ------------|-----------------
!>  0.80        |   1.281551565545
!>  0.90        |   1.644853626951
!>  0.95        |   1.959963984540
!>  0.98        |   2.326347874041
!>  0.99        |   2.575829303549
!>  0.995       |   2.807033768344
!>  0.998       |   3.090232306168
!>  0.999       |   3.290526731492
!>  0.9999      |   3.890591886413
!>  0.99999     |   4.417173413469
!>  0.999999    |   4.891638475699
!>  0.9999999   |   5.326723886384
!>  0.99999999  |   5.730728868236
!>  0.999999999 |   6.109410204869
!>
!>  For small \f$p\f$, the quantile function has the useful asymptotic expansion \f$\Phi^{-1}(p) = -\sqrt{\ln{\frac{1}{p^{2}}} - \ln\ln{\frac{1}{p^{2}}} - \ln(2\pi)} + \mathcal{o}(1)\f$.
!>
!>  **Random Number Generation**<br>
!>
!>  The current implementations of the RNG generic interfaces of this module use
!>  the Box-Muller trigonometric and rejection methods for Normal random number generation.<br>
!>
!>  \note
!>  The `real32` Standard Normal random numbers generated by the Box-Muller algorithm are known to be limited to the range \f$[-6.66, +6.66]\f$.<br>
!>  The `real64` Standard Normal random numbers generated by the Box-Muller algorithm are known to be limited to the range \f$[-9.419, +9.419]\f$.<br>
!>  **Entropy**<br>
!>
!>  The **entropy** of the Normal distribution is defined by the following equation,
!>  \f{equation}{
!>      \large
!>      \mathcal{H}(\sigma^2) = \frac{1}{2} \log(2\pi\sigma^2) + \frac{1}{2}
!>  \f}
!>  **Fisher Information**<br>
!>
!>  The Fisher information for the Normal distribution is defined by the following equation,
!>  \f{equation}{
!>      \large
!>      \mathcal{I}(\mu,\sigma) =
!>      \begin{pmatrix}
!>          \frac{1}{\sigma^2} & 0 \\
!>          0 & \frac{2}{\sigma^2} \\
!>      \end{pmatrix}
!>  \f}
!>
!>  **Kullback-Leibler Divergence (KLD)**<br>
!>
!>  The Kullback-Leibler Divergence, also known as the **relative entropy**, of a univariate Normal
!>  distribution \f$Q\f$ from a reference univariate Normal distribution \f$P\f$ is defined as,
!>  \f{equation}{
!>      \large
!>      D_{KL}(P \parallel Q) =
!>      \frac{(\mu_P - \mu_Q)^2}{2\sigma_Q^2} +
!>      \frac{1}{2}\bigg( \frac{\sigma_P}{\sigma_Q} - \ln\big(\frac{\sigma_P}{\sigma_Q}\big) - 1 \bigg) ~,
!>  \f}
!>
!>  \see
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{ziggurat_vs_boxbasic, The runtime performance of [setNormRand](@ref pm_distNorm::setNormRand) and [setNormRandBox](@ref pm_distNorm::setNormRandBox).}
!>  \include{lineno} benchmark/pm_distNorm/ziggurat_vs_boxbasic/main.F90
!>  \compilefb{ziggurat_vs_boxbasic}
!>  \postprocb{ziggurat_vs_boxbasic}
!>  \include{lineno} benchmark/pm_distNorm/ziggurat_vs_boxbasic/main.py
!>  \visb{ziggurat_vs_boxbasic}
!>  \image html benchmark/pm_distNorm/ziggurat_vs_boxbasic/benchmark.ziggurat_vs_boxbasic.runtime.png width=1000
!>  \image html benchmark/pm_distNorm/ziggurat_vs_boxbasic/benchmark.ziggurat_vs_boxbasic.runtime.ratio.png width=1000
!>  \moralb{ziggurat_vs_boxbasic}
!>      -#  The benchmark procedures named `setNormRandZiggurat` and `setNormRandBoxBasic` call the generic interfaces [setNormRand](@ref pm_distNorm::setNormRand)
!>          (ziggurat method) and [setNormRandBox](@ref pm_distNorm::setNormRandBox) (basic trigonometric Box-Muller method) respectively.<br>
!>      -#  The benchmark procedure named `setNormRandZigX256S` calls the generic interface [setNormRand](@ref pm_distNorm::setNormRand)
!>          with the [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) RNG instead of the default intrinsic Fortran RNG.<br>
!>      -#  While [setNormRandBox](@ref pm_distNorm::setNormRandBox) are pure subroutines,
!>          [setNormRand](@ref pm_distNorm::setNormRand) is `impure` when the default intrinsic Fortran RNG is used.<br>
!>      -#  The benchmark results confirm the widespread community observation that the Ziggurat method for Normal RNG can slightly faster than the Box-Muller method.<br>
!>          The difference, however, appears to be marginal, compiler and hardware dependent.<br>
!>
!>  \see
!>  [pm_distLogNorm](@ref pm_distLogNorm)<br>
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  [pm_distNormShell](@ref pm_distNormShell)<br>
!>
!>  \test
!>  [test_pm_distNorm](@ref test_pm_distNorm)
!>
!>  \todo
!>  \pvlow
!>  The performance of the current implementation of the Ziggurat RNG algorithm can be slightly improved for single-precision Normal RNG
!>  in the procedures that rely on [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) for uniform RNG.<br>
!>  This requires a customized implementation of the Ziggurat method for this `real` kind.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distNorm

    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_kind, only: SK, IK, LK, RKH, RKB

    implicit none

    character(*, SK), parameter         :: MODULE_NAME = "@pm_distNorm"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Univariate Normal
    !>  as defined in the description of [pm_distNorm](@ref pm_distNorm).
    !>
    !>  \details
    !>  See the documentation of [pm_distNorm](@ref pm_distNorm) for the definition of the Univariate Normal distribution.
    !>
    !>  \interface{distNorm_type}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: distNorm_type
    !>      type(distNorm_type) :: distNorm
    !>
    !>      distNorm = distNorm_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distNorm_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distNorm_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the natural logarithm of probability density function (PDF) of the univariate Normal distribution.
    !>
    !>  \details
    !>  See the documentation of [pm_distNorm](@ref pm_distNorm) for further details.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `logPDF`, representing the point(s) at which the PDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as the output `logPDF` representing the mean of the distribution.<br>
    !>                              (**optional**, default = `0.`)
    !>  \param[in]  sigma       :   The input scalar of the same type and kind as the output `logPDF`
    !>                              representing the inverse of the standard deviation of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as the input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              representing the natural logarithm of the PDF of the Normal distribution at the specified `x`.
    !>
    !>  \interface{getNormLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormLogPDF
    !>
    !>      logPDF = getNormLogPDF(x, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding procedure argument.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  See [setNormLogPDF](@ref pm_distNorm::setNormLogPDF) for a more
    !>  performant but less flexible interface of the same functionality.<br>
    !>
    !>  \see
    !>  [setNormLogPDF](@ref pm_distNorm::setNormLogPDF)<br>
    !>  [getLogNormLogPDF](@ref pm_distLogNorm::getLogNormLogPDF)<br>
    !>  [setLogNormLogPDF](@ref pm_distLogNorm::setLogNormLogPDF)<br>
    !>
    !>  \example{getNormLogPDF}
    !>  \include{lineno} example/pm_distNorm/getNormLogPDF/main.F90
    !>  \compilef{getNormLogPDF}
    !>  \output{getNormLogPDF}
    !>  \include{lineno} example/pm_distNorm/getNormLogPDF/main.out.F90
    !>  \postproc{getNormLogPDF}
    !>  \include{lineno} example/pm_distNorm/getNormLogPDF/main.py
    !>  \vis{getNormLogPDF}
    !>  \image html pm_distNorm/getNormLogPDF/getNormLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{getNormLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNormLogPDF_RK5(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNormLogPDF_RK4(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNormLogPDF_RK3(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNormLogPDF_RK2(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNormLogPDF_RK1(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormLogPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the natural logarithm of probability density function (PDF) of the univariate Normal distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar or array of the same shape as the input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              representing the natural logarithm of the PDF of the Normal distribution at `x`.<br>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `logPDF`, representing the point(s) at which the PDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type and
    !>                              kind as the output `logPDF` representing the mean of the distribution.<br>
    !>                              (**optional**, default = `0.`)
    !>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `logPDF` representing the inverse of the standard deviation of the distribution.<br>
    !>                              (**optional**, default = `1.`, must be present <b>if and only if</b> `logInvSigma` is also present.)
    !>  \param[in]  logInvSigma :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `logPDF` representing the natural logarithm of the inverse of the standard
    !>                              deviation of the distribution.<br>
    !>                              (**optional**, default = `0`, must be present <b>if and only if</b> `invSigma` is also present.)
    !>
    !>  \interface{setNormLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: setNormLogPDF
    !>
    !>      call setNormLogPDF(logPDF, x)
    !>      call setNormLogPDF(logPDF, x, mu)
    !>      call setNormLogPDF(logPDF, x, invSigma, logInvSigma)
    !>      call setNormLogPDF(logPDF, x, mu, invSigma, logInvSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < invSigma` must hold for the corresponding procedure argument.<br>
    !>  The condition `log(invSigma) == logInvSigma` must hold within a small range of computer precision
    !>  for the corresponding procedure arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  See [getNormLogPDF](@ref pm_distNorm::getNormLogPDF) for a less
    !>  performant but more flexible interface of the same functionality.<br>
    !>
    !>  \see
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [getLogNormLogPDF](@ref pm_distLogNorm::getLogNormLogPDF)<br>
    !>  [setLogNormLogPDF](@ref pm_distLogNorm::setLogNormLogPDF)<br>
    !>
    !>  \example{setNormLogPDF}
    !>  \include{lineno} example/pm_distNorm/setNormLogPDF/main.F90
    !>  \compilef{setNormLogPDF}
    !>  \output{setNormLogPDF}
    !>  \include{lineno} example/pm_distNorm/setNormLogPDF/main.out.F90
    !>  \postproc{setNormLogPDF}
    !>  \include{lineno} example/pm_distNorm/setNormLogPDF/main.py
    !>  \vis{setNormLogPDF}
    !>  \image html pm_distNorm/setNormLogPDF/setNormLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \todo
    !>  \pmed
    !>  A performant vectorized `logPDF(:)` version of the subroutines under this generic interface could be added in the future.
    !>
    !>  \final{setNormLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNormLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormLogPDFDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormLogPDFDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormLogPDFDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormLogPDFDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormLogPDFDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormLogPDFMD_RK5(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormLogPDFMD_RK4(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormLogPDFMD_RK3(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormLogPDFMD_RK2(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormLogPDFMD_RK1(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormLogPDFDS_RK5(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormLogPDFDS_RK4(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormLogPDFDS_RK3(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormLogPDFDS_RK2(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormLogPDFDS_RK1(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFDS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormLogPDFMS_RK5(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormLogPDFMS_RK4(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormLogPDFMS_RK3(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormLogPDFMS_RK2(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormLogPDFMS_RK1(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormLogPDFMS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setNormLogPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the univariate Normal distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `cdf`, representing the point(s) at which the CDF must be computed.
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as the output `cdf` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  sigma       :   The input scalar of the same type and kind as the output `cdf`
    !>                              representing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as the input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the CDF of the distribution at the specified input `x`.
    !>
    !>  \interface{getNormCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormCDF
    !>
    !>      cdf = getNormCDF(x, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding procedure argument.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setNormCDF](@ref pm_distNorm::setNormCDF)<br>
    !>
    !>  \example{getNormCDF}
    !>  \include{lineno} example/pm_distNorm/getNormCDF/main.F90
    !>  \compilef{getNormCDF}
    !>  \output{getNormCDF}
    !>  \include{lineno} example/pm_distNorm/getNormCDF/main.out.F90
    !>  \postproc{getNormCDF}
    !>  \include{lineno} example/pm_distNorm/getNormCDF/main.py
    !>  \vis{getNormCDF}
    !>  \image html pm_distNorm/getNormCDF/getNormCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{getNormCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNormCDF_RK5(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNormCDF_RK4(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNormCDF_RK3(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNormCDF_RK2(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNormCDF_RK1(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the univariate Normal distribution.
    !>
    !>  \param[out] cdf         :   The output scalar or array of the same shape as the input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the CDF of the distribution at the specified input `x`.
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `cdf`, representing the point(s) at which the CDF must be computed.
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as the output `cdf` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  invSigma    :   The input scalar of the same type and kind as the output `quantile`,
    !>                              representing the **inverse** of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>
    !>  \interface{setNormCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: setNormCDF
    !>
    !>      call setNormCDF(cdf, x)
    !>      call setNormCDF(cdf, x, mu)
    !>      call setNormCDF(cdf, x, mu, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `invSigma > 0.` must hold for the corresponding procedure argument.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getNormCDF](@ref pm_distNorm::getNormCDF)<br>
    !>
    !>  \example{setNormCDF}
    !>  \include{lineno} example/pm_distNorm/setNormCDF/main.F90
    !>  \compilef{setNormCDF}
    !>  \output{setNormCDF}
    !>  \include{lineno} example/pm_distNorm/setNormCDF/main.out.F90
    !>  \postproc{setNormCDF}
    !>  \include{lineno} example/pm_distNorm/setNormCDF/main.py
    !>  \vis{setNormCDF}
    !>  \image html pm_distNorm/setNormCDF/setNormCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{setNormCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNormCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormCDFDD_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormCDFDD_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormCDFDD_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormCDFDD_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(out)                           :: cdf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormCDFDD_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormCDFMD_RK5(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormCDFMD_RK4(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormCDFMD_RK3(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormCDFMD_RK2(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormCDFMD_RK1(cdf, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormCDFMS_RK5(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormCDFMS_RK4(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormCDFMS_RK3(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormCDFMS_RK2(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormCDFMS_RK1(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormCDFMS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Normal Quantile corresponding to the input CDF of the univariate Normal distribution.
    !>
    !>  \param[in]  cdf         :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `quantile`, representing the point(s) (probabilities) at which the quantile must be computed.
    !>                              The `cdf` must be in the range \f$[0, 1]\f$.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as the output `quantile` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  sigma       :   The input scalar of the same type and kind as the output `quantile`,
    !>                              representing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `quantile`              :   The output scalar or array of the same shape as the input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the quantile of the distribution at the specified input `cdf`.
    !>
    !>  \interface{getNormQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormQuan
    !>
    !>      quantile = getNormQuan(cdf, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding procedure arguments.<br>
    !>  The condition `0. <= cdf .and. cdf <= 1.` must hold for the corresponding procedure arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setNormQuan](@ref pm_distNorm::setNormQuan)<br>
    !>
    !>  \example{getNormQuan}
    !>  \include{lineno} example/pm_distNorm/getNormQuan/main.F90
    !>  \compilef{getNormQuan}
    !>  \output{getNormQuan}
    !>  \include{lineno} example/pm_distNorm/getNormQuan/main.out.F90
    !>  \postproc{getNormQuan}
    !>  \include{lineno} example/pm_distNorm/getNormQuan/main.py
    !>  \vis{getNormQuan}
    !>  \image html pm_distNorm/getNormQuan/getNormQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{getNormQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNormQuan_RK5(cdf, mu, sigma) result(quantile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormQuan_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: quantile
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNormQuan_RK4(cdf, mu, sigma) result(quantile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormQuan_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: quantile
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNormQuan_RK3(cdf, mu, sigma) result(quantile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormQuan_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: quantile
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNormQuan_RK2(cdf, mu, sigma) result(quantile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormQuan_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: quantile
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNormQuan_RK1(cdf, mu, sigma) result(quantile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormQuan_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: quantile
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormQuan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the quantile of the univariate Normal distribution at the specified input CDF.
    !>
    !>  \param[out] quantile    :   The output scalar or array of the same shape as the input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the quantile of the distribution at the specified input `cdf`.<br>
    !>  \param[in]  cdf         :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `quantile`, representing the point(s) (probabilities) at which the quantile must be computed.<br>
    !>                              The `cdf` must be in the range \f$[0, 1]\f$.<br>
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as the output `quantile` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  sigma       :   The input scalar of the same type and kind as the output `quantile`,
    !>                              representing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>
    !>  \interface{setNormQuan}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: setNormQuan
    !>
    !>      call setNormQuan(quantile, cdf)
    !>      call setNormQuan(quantile, cdf, mu)
    !>      call setNormQuan(quantile, cdf, mu, sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding procedure argument.<br>
    !>  The condition `0. <= cdf .and. cdf <= 1.` must hold for the corresponding procedure arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The procedures of this generic interface are merely wrappers around the procedures of [setErfInv](@ref pm_mathErf::setErfInv).<br>
    !>  See the documentation of [setErfInv](@ref pm_mathErf::setErfInv) about the accuracy of the computed quantile.<br>
    !>
    !>  \see
    !>  [getNormQuan](@ref pm_distNorm::getNormQuan)<br>
    !>
    !>  \example{setNormQuan}
    !>  \include{lineno} example/pm_distNorm/setNormQuan/main.F90
    !>  \compilef{setNormQuan}
    !>  \output{setNormQuan}
    !>  \include{lineno} example/pm_distNorm/setNormQuan/main.out.F90
    !>  \postproc{setNormQuan}
    !>  \include{lineno} example/pm_distNorm/setNormQuan/main.py
    !>  \vis{setNormQuan}
    !>  \image html pm_distNorm/setNormQuan/setNormQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{setNormQuan}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNormQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormQuanDD_RK5(quantile, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormQuanDD_RK4(quantile, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormQuanDD_RK3(quantile, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormQuanDD_RK2(quantile, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(out)                           :: quantile
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormQuanDD_RK1(quantile, cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormQuanMD_RK5(quantile, cdf, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormQuanMD_RK4(quantile, cdf, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormQuanMD_RK3(quantile, cdf, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormQuanMD_RK2(quantile, cdf, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormQuanMD_RK1(quantile, cdf, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormQuanMS_RK5(quantile, cdf, mu, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormQuanMS_RK4(quantile, cdf, mu, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormQuanMS_RK3(quantile, cdf, mu, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormQuanMS_RK2(quantile, cdf, mu, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormQuanMS_RK1(quantile, cdf, mu, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormQuanMS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: quantile
        real(RKG)           , intent(in)                            :: cdf
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setNormQuan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `integer` of default kind \IK containing the output of Fortran intrinsic `precision()`
    !>  for the `real` kind used to generate the constant array [ZIG_RKB](@ref pm_distNorm::ZIG_RKB).<br>
    !>
    !>  \see
    !>  [setZig](@ref pm_ziggurat::getZig)
    !>  [ZIG_RKB](@ref pm_distNorm::ZIG_RKB)
    !>  [setNormRand](@ref pm_distNorm::setNormRand)
    !>  [ZIG_PRECISION](@ref pm_distNorm::ZIG_PRECISION)
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{ZIG_PRECISION}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    integer(IK) , parameter :: ZIG_PRECISION = 33_IK

    real(RKB), parameter, private :: ZIGSET1(2, 0:128) = reshape([ +3.910757959524915869549621434510571640_RKB, +0.000000000000000000000000000000000000E+0_RKB &
                                                                , +3.654152885361008771645429720399515670_RKB, +0.502781343070952005087938400241541891E-3_RKB &
                                                                , +3.449278298561431270627228213833611250_RKB, +0.104086943106322506013635915064907845E-2_RKB &
                                                                , +3.320244733839825517532232984442230700_RKB, +0.161091799459473453538408742342021986E-2_RKB &
                                                                , +3.224575052047801587144019828764775500_RKB, +0.220312016549958827311372671027399889E-2_RKB &
                                                                , +3.147889289518000685451855194084413580_RKB, +0.281289233937616172611348633870502865E-2_RKB &
                                                                , +3.083526132002143251877768947617198840_RKB, +0.343751917929162224416229759239016324E-2_RKB &
                                                                , +3.027837791769593524571714584215607150_RKB, +0.407518400039001087464038448139061472E-2_RKB &
                                                                , +2.978603279881843165536974212294808030_RKB, +0.472457682607575925888895096351016103E-2_RKB &
                                                                , +2.934366867208887589959928979567677640_RKB, +0.538470372266379748395288853907771468E-2_RKB &
                                                                , +2.894121053613412181388100356210103800_RKB, +0.605478221942164063991037466747021604E-2_RKB &
                                                                , +2.857138730873224588561645268053226920_RKB, +0.673417886624137640800598237565554838E-2_RKB &
                                                                , +2.822877396826442907534115515646593530_RKB, +0.742236950888280249329920345563351449E-2_RKB &
                                                                , +2.790921174001927318997779045468293070_RKB, +0.811891273868296054077142626445965006E-2_RKB &
                                                                , +2.760944005279986201244382392492695240_RKB, +0.882343143223528339445647449425005264E-2_RKB &
                                                                , +2.732685359044011420043182513046678140_RKB, +0.953559949343287212693559957627379651E-2_RKB &
                                                                , +2.705933656123062221333700225998574460_RKB, +0.102551320707227866483396053794922990E-1_RKB &
                                                                , +2.680514643285745101098374611431323380_RKB, +0.109817781711791129427823736030159900E-1_RKB &
                                                                , +2.656283037576743296802124304559512000_RKB, +0.117153149732408013682731286377880085E-1_RKB &
                                                                , +2.633116393631582759976309292516893140_RKB, +0.124555433719091587849084452906815838E-1_RKB &
                                                                , +2.610910518488823671930263694478955450_RKB, +0.132022844366303809959357880408105012E-1_RKB &
                                                                , +2.589575986708286649808805574507900700_RKB, +0.139553765573260408177619288977468218E-1_RKB &
                                                                , +2.569035452681843781314262921528081910_RKB, +0.147146731176153334234423656007745630E-1_RKB &
                                                                , +2.549221550324783104422671371241746650_RKB, +0.154800405777112059806768992081788185E-1_RKB &
                                                                , +2.530075232159854187716539084126876370_RKB, +0.162513568797661774124156785927309438E-1_RKB &
                                                                , +2.511544441626694343254607344580084390_RKB, +0.170285101099635321775373808069033731E-1_RKB &
                                                                , +2.493583041271046768170053296506813290_RKB, +0.178113973671931254109163549358583312E-1_RKB &
                                                                , +2.476149939670523163756216268162166500_RKB, +0.185999237995485118944697972790344771E-1_RKB &
                                                                , +2.459208374334705035673859596487019180_RKB, +0.193940017783549904963609637793049119E-1_RKB &
                                                                , +2.442725318200364223794234919225453500_RKB, +0.201935501858165129428739457647187713E-1_RKB &
                                                                , +2.426670984937146719863529851633726960_RKB, +0.209984937972267218821328887239805004E-1_RKB &
                                                                , +2.411018413901119491690349211725038630_RKB, +0.218087627424279111659346610629834973E-1_RKB &
                                                                , +2.395743119781927356168686681412574760_RKB, +0.226242920341075527766680294896069533E-1_RKB &
                                                                , +2.380822795172085556506619691319380400_RKB, +0.234450211528013514383037474818822027E-1_RKB &
                                                                , +2.366237056717290911362148128140186000_RKB, +0.242708936802748765008871090362921378E-1_RKB &
                                                                , +2.351967227379144761902530751453964640_RKB, +0.251018569743934986305835005740321359E-1_RKB &
                                                                , +2.337996148796528635433480327093713940_RKB, +0.259378618797451777000075273655017443E-1_RKB &
                                                                , +2.324308018871132508266119157050790040_RKB, +0.267788624692146997471273278024758566E-1_RKB &
                                                                , +2.310888250601371758550614355860895420_RKB, +0.276248158124683715174543699462393572E-1_RKB &
                                                                , +2.297723348902863520079790814230529130_RKB, +0.284756817679310236840265332909329685E-1_RKB &
                                                                , +2.284800802724492127387834486937223350_RKB, +0.293314227953502809501556980682789090E-1_RKB &
                                                                , +2.272108990228381861937683717373910900_RKB, +0.301920037864680654091842729977054983E-1_RKB &
                                                                , +2.259637095173787624597566531167311770_RKB, +0.310573919116731481152577086508924965E-1_RKB &
                                                                , +2.247375032947389262297952392515659590_RKB, +0.319275564808046037696524061552362067E-1_RKB &
                                                                , +2.235313384929921110748362199673550240_RKB, +0.328024688165248443463852626030938600E-1_RKB &
                                                                , +2.223443340092510611365346409716922540_RKB, +0.336821021388909593274737975695380343E-1_RKB &
                                                                , +2.211756642884160997470500709266731210_RKB, +0.345664314599311541791934500133628777E-1_RKB &
                                                                , +2.200245546611276427712165173291835590_RKB, +0.354554334871846263278378901217329216E-1_RKB &
                                                                , +2.188902771626360742839576505604305880_RKB, +0.363490865352926868627480678085691155E-1_RKB &
                                                                , +2.177721467740293002579164079152503390_RKB, +0.372473704448399433976210004282028638E-1_RKB &
                                                                , +2.166695180354308542353137142121026540_RKB, +0.381502665077398667348211240367435477E-1_RKB &
                                                                , +2.155817819876737469119503677197502650_RKB, +0.390577573985415074098963321878992358E-1_RKB &
                                                                , +2.145083634047888982767999729572053870_RKB, +0.399698271111055203479683284323823695E-1_RKB &
                                                                , +2.134487182846016909178836604796243250_RKB, +0.408864609001596660576863760104051322E-1_RKB &
                                                                , +2.124023315689523545420714787478384860_RKB, +0.418076452272979751244992871255722388E-1_RKB &
                                                                , +2.113687150686653177781935198589005610_RKB, +0.427333677110349512002228670090424488E-1_RKB &
                                                                , +2.103474055714877305933714304045440450_RKB, +0.436636170805675209921984193602696278E-1_RKB &
                                                                , +2.093379631138791930166361585961667970_RKB, +0.445983831329337385181389634688869380E-1_RKB &
                                                                , +2.083399693998304613670792088175897210_RKB, +0.455376566932892082093887264381438529E-1_RKB &
                                                                , +2.073530263518743034646393248416971800_RKB, +0.464814295780503946152953311870656232E-1_RKB &
                                                                , +2.063767547811732114341853749078149630_RKB, +0.474296945606789346776298063351001244E-1_RKB &
                                                                , +2.054107931650652130219475666313816510_RKB, +0.483824453399031873773625284851490337E-1_RKB &
                                                                , +2.044547965217531455282628792749082960_RKB, +0.493396765101929088581344120097015200E-1_RKB &
                                                                , +2.035084353729618971413871948694988270_RKB, +0.503013835343204408492863489427518069E-1_RKB &
                                                                , +2.025713947863854245252399025955272060_RKB, +0.512675627178574145127079248497036520E-1_RKB &
                                                                , +2.016433734906204123873988577023278720_RKB, +0.522382111854699318500836298767666937E-1_RKB &
                                                                , +2.007240830560528758913239738198299090_RKB, +0.532133268588876922360496499106520185E-1_RKB &
                                                                , +1.998132471358419680392214769162858870_RKB, +0.541929084364337554633405727268610055E-1_RKB &
                                                                , +1.989106007617438123201023091217180910_RKB, +0.551769553740117250736639942272910314E-1_RKB &
                                                                , +1.980158896900476605540416694815765850_RKB, +0.561654678674562273820736142229565223E-1_RKB &
                                                                , +1.971288697933659294606356474434686550_RKB, +0.571584468361607664497222499802473837E-1_RKB &
                                                                , +1.962493064944363052826024381256974850_RKB, +0.581558939079044529261453899150298019E-1_RKB &
                                                                , +1.953769742384646776692572811817963060_RKB, +0.591578114048058224055549195846209505E-1_RKB &
                                                                , +1.945116560008678301234686281875394240_RKB, +0.601642023303380532839147151743642514E-1_RKB &
                                                                , +1.936531428275694700380702197135648770_RKB, +0.611750703573454323914419246868044098E-1_RKB &
                                                                , +1.928012334052665710328808957510011910_RKB, +0.621904198170059582295143354352838012E-1_RKB &
                                                                , +1.919557336593188113064282063465457130_RKB, +0.632102556886895688287062083850491706E-1_RKB &
                                                                , +1.911164563771253338349016119398665420_RKB, +0.642345835906656803805852580678544053E-1_RKB &
                                                                , +1.902832208550429269452784399230803680_RKB, +0.652634097716175649174146849378032930E-1_RKB &
                                                                , +1.894558525670704732040885690575958130_RKB, +0.662967411029246168483446165075442866E-1_RKB &
                                                                , +1.886341828536782820037856553309563160_RKB, +0.673345850716767915010299444811729445E-1_RKB &
                                                                , +1.878180486292995844684467996043645020_RKB, +0.683769497743884728203642791575539760E-1_RKB &
                                                                , +1.870072921071266778496337776290841830_RKB, +0.694238439113817677925134633773045847E-1_RKB &
                                                                , +1.862017605399674118665348784397922740_RKB, +0.704752767818117550133248257276187773E-1_RKB &
                                                                , +1.854013059760201906750898441184461490_RKB, +0.715312582793085547211649297744966411E-1_RKB &
                                                                , +1.846057850285185505570580952991448920_RKB, +0.725917988882132560598516526262197654E-1_RKB &
                                                                , +1.838150586582806633764931304898561700_RKB, +0.736569096803867509488343563799016629E-1_RKB &
                                                                , +1.830289919682756933756055676984081970_RKB, +0.747266023125723976798301015314803279E-1_RKB &
                                                                , +1.822474540093885838871898571617810260_RKB, +0.758008890242951847298446370142010715E-1_RKB &
                                                                , +1.814703175966282671680772343522476780_RKB, +0.768797826362816984781690401687711081E-1_RKB &
                                                                , +1.806974591350820938703420428388107530_RKB, +0.779632965493867285846698800079747914E-1_RKB &
                                                                , +1.799287584549720199341728680984235510_RKB, +0.790514447440137817496961412302243116E-1_RKB &
                                                                , +1.791640986552162594624382322132709310_RKB, +0.801442417800181275436054931177368005E-1_RKB &
                                                                , +1.784033659549441512971864077807943680_RKB, +0.812417027970822772687529532474258376E-1_RKB &
                                                                , +1.776464495524522868996124119276484200_RKB, +0.823438435155550059817881401744955979E-1_RKB &
                                                                , +1.768932414911268589029665159840958690_RKB, +0.834506802377461758021915029242417507E-1_RKB &
                                                                , +1.761436365318910280539794251593908460_RKB, +0.845622298496707118376316534786168010E-1_RKB &
                                                                , +1.753975320317671535176286495650271100_RKB, +0.856785098232361263350014624915617177E-1_RKB &
                                                                , +1.746548278281722412853610819915111290_RKB, +0.867995382188689874348253633476266673E-1_RKB &
                                                                , +1.739154261285911657262420055628115450_RKB, +0.879253336885766911833903717911534395E-1_RKB &
                                                                , +1.731792314052963154137933505808925000_RKB, +0.890559154794418239043567001488716945E-1_RKB &
                                                                , +1.724461502948044912052862397775662270_RKB, +0.901913034375473009991491331657262247E-1_RKB &
                                                                , +1.717160915017823089741659119609676780_RKB, +0.913315180123313418058911041645338441E-1_RKB &
                                                                , +1.709889657071301820241745482016102880_RKB, +0.924765802613721921296195811766027943E-1_RKB &
                                                                , +1.702646854799923151653900407051805770_RKB, +0.936265118556033400788781021199070560E-1_RKB &
                                                                , +1.695431651934561568299588089192625460_RKB, +0.947813350849607903374772538857282822E-1_RKB &
                                                                , +1.688243209437195389093695018261921410_RKB, +0.959410728644647702377950152378611650E-1_RKB &
                                                                , +1.681080704725173871909617303396002540_RKB, +0.971057487407390411201725488177809241E-1_RKB &
                                                                , +1.673943330926124999231911729230294350_RKB, +0.982753868989717834844303967803275557E-1_RKB &
                                                                , +1.666830296161665512280717200364773540_RKB, +0.994500121703228172933039137854161328E-1_RKB &
                                                                , +1.659740822858182552384789747414556960_RKB, +0.100629650039782712328127609576129162E+0_RKB &
                                                                , +1.652674147083055944977076666279072050_RKB, +0.101814326654490140522320621972065562E+0_RKB &
                                                                , +1.645629517904782346099997435068184430_RKB, +0.103004068832514625466733060901486566E+0_RKB &
                                                                , +1.638606196775547730191205373889528150_RKB, +0.104198904072112656528428960098709859E+0_RKB &
                                                                , +1.631603456934873546471532683464361160_RKB, +0.105398860561465958979782511482284346E+0_RKB &
                                                                , +1.624620582833034778354059455865975410_RKB, +0.106603967188911549935198869215607602E+0_RKB &
                                                                , +1.617656869573015532637880619762087450_RKB, +0.107814253553674065499513686193278908E+0_RKB &
                                                                , +1.610711622369830051160545249046193620_RKB, +0.109029749977111720145296612324748632E+0_RKB &
                                                                , +1.603784156026094530393066805218239570_RKB, +0.110250487514488157810061008578271138E+0_RKB &
                                                                , +1.596873794422788175563087978818873180_RKB, +0.111476497967283378462230486935475436E+0_RKB &
                                                                , +1.589979870024190797461111747861209560_RKB, +0.112707813896057876934855393386708534E+0_RKB &
                                                                , +1.583101723396029247521107058058590910_RKB, +0.113944468633885115777321981289805190E+0_RKB &
                                                                , +1.576238702735906320876658002975291310_RKB, +0.115186496300368473953376646626857215E+0_RKB &
                                                                , +1.569390163415123656042832370955675890_RKB, +0.116433931816259871778263682344311485E+0_RKB &
                                                                , +1.562555467531044820930627254229273630_RKB, +0.117686810918698373049736317381864947E+0_RKB &
                                                                , +1.555733983469176375991038619565492320_RKB, +0.118945170177088211568318678970235501E+0_RKB &
                                                                , +1.548925085474173406375274133330557850_RKB, +0.120209047009636885032310150781360772E+0_RKB &
                                                                , +1.542128153229001959318075945656599760_RKB, +0.121478479700575208713903649777558606E+0_RKB &
                                                                , +1.535342571441514138082707130216935300_RKB, +0.122753507418082528688239173356532507E+0_RKB ], shape = [2, 129])
    real(RKB), parameter, private :: ZIGSET2(2, 1:128) = reshape([ +1.528567729437712402628199461442167180_RKB, +0.124034170232941664267186768442998377E+0_RKB &
                                                                , +1.521803020760998008679228498636570410_RKB, +0.125320509137949586535569901372920916E+0_RKB &
                                                                , +1.515047842776714566947194255734634790_RKB, +0.126612566068111349659514010876108850E+0_RKB &
                                                                , +1.508301596281311496171425183785073880_RKB, +0.127910383921646379432792471520964754E+0_RKB &
                                                                , +1.501563685115463738688173273682448190_RKB, +0.129214006581837895215365859576747957E+0_RKB &
                                                                , +1.494833515780493555035890159317754310_RKB, +0.130523478939758003270817270840923056E+0_RKB &
                                                                , +1.488110497057447553013103271226758760_RKB, +0.131838846917902858239534927650542548E+0_RKB &
                                                                , +1.481394039628187363902263028944763970_RKB, +0.133160157494774252287691074366070633E+0_RKB &
                                                                , +1.474683555697855570628865196528769660_RKB, +0.134487458730446066069814161355648980E+0_RKB &
                                                                , +1.467978458618079624962505551018572310_RKB, +0.135820799793156210330787895235534157E+0_RKB &
                                                                , +1.461278162510275558141634572902431380_RKB, +0.137160230986967010673991854114446593E+0_RKB &
                                                                , +1.454582081888410275202473614895370850_RKB, +0.138505803780539450342898439610712141E+0_RKB &
                                                                , +1.447889631280576100338960603208179050_RKB, +0.139857570837069297157912765350925466E+0_RKB &
                                                                , +1.441200224848723969896483905460819290_RKB, +0.141215586045435912189035654183361969E+0_RKB &
                                                                , +1.434513276005892200373197821347744230_RKB, +0.142579904552617481390538701477587025E+0_RKB &
                                                                , +1.427828197030256028046923580264626420_RKB, +0.143950582797429540314021584008091148E+0_RKB &
                                                                , +1.421144398675309048678328840045534000_RKB, +0.145327678545646990255525094264574683E+0_RKB &
                                                                , +1.414461289775471190729091975769526810_RKB, +0.146711250926573347052911625905291106E+0_RKB &
                                                                , +1.407778276846398829890729544114189010_RKB, +0.148101360471124737782233265331001140E+0_RKB &
                                                                , +1.401094763679250977372468027391928780_RKB, +0.149498069151500183758164830605819877E+0_RKB &
                                                                , +1.394410150928141013910818209046742880_RKB, +0.150901440422514000012107606564276651E+0_RKB &
                                                                , +1.387723835689976042816457610158740000_RKB, +0.152311539264670722976597120096621794E+0_RKB &
                                                                , +1.381035211075855426557023886048995030_RKB, +0.153728432229067872471998597962422358E+0_RKB &
                                                                , +1.374343665773166259809330405021183000_RKB, +0.155152187484217086331763054685730862E+0_RKB &
                                                                , +1.367648583597476202662733878746826080_RKB, +0.156582874864879763414378129073938725E+0_RKB &
                                                                , +1.360949343033283011396528844287145660_RKB, +0.158020565923019343095774015163019653E+0_RKB &
                                                                , +1.354245316762634995007843766072061330_RKB, +0.159465333980978769091659649295558937E+0_RKB &
                                                                , +1.347535871180587198226303191916842980_RKB, +0.160917254186998568094181273666041719E+0_RKB &
                                                                , +1.340820365896404038797740514363138740_RKB, +0.162376403573198357993539454690214112E+0_RKB &
                                                                , +1.334098153219360045667959181721837720_RKB, +0.163842861116152528813630294472766962E+0_RKB &
                                                                , +1.327368577627925853644432080219569470_RKB, +0.165316707800199358376684187822216740E+0_RKB &
                                                                , +1.320630975221056264316363281880977050_RKB, +0.166798026683631985045755192821871952E+0_RKB &
                                                                , +1.313884673150220489854736808952840970_RKB, +0.168286902967929517543098016930008494E+0_RKB &
                                                                , +1.307128989030731110178131195502915220_RKB, +0.169783424070197178160245257400787177E+0_RKB &
                                                                , +1.300363230330837190351929439655934200_RKB, +0.171287679698995818105358886788417999E+0_RKB &
                                                                , +1.293586693736947753945601406359300470_RKB, +0.172799761933753486486498659023686530E+0_RKB &
                                                                , +1.286798664493243646316660896517199420_RKB, +0.174319765307965059246455246546117697E+0_RKB &
                                                                , +1.279998415713817924797308907917922830_RKB, +0.175847786896400331371952020571991623E+0_RKB &
                                                                , +1.273185207665356364764060323720764170_RKB, +0.177383926406556544373322976411406151E+0_RKB &
                                                                , +1.266358287018229453775313848849852450_RKB, +0.178928286274608171283092549748417432E+0_RKB &
                                                                , +1.259516886063714228190559472001646690_RKB, +0.180480971766125034830178787412733260E+0_RKB &
                                                                , +1.252660221894897227354473865616388090_RKB, +0.182042091081849625639024292847696263E+0_RKB &
                                                                , +1.245787495548627294596662055366290760_RKB, +0.183611755468845965529740751119222664E+0_RKB &
                                                                , +1.238897891105687374493975672708951630_RKB, +0.185190079337355691907687938603429913E+0_RKB &
                                                                , +1.231990574746136091354596183983300960_RKB, +0.186777180383722406901836668848789685E+0_RKB &
                                                                , +1.225064693756530787096859391988721310_RKB, +0.188373179719772944131913560975705614E+0_RKB &
                                                                , +1.218119375485481656036492926236470550_RKB, +0.189978202009074284878057268848586717E+0_RKB &
                                                                , +1.211153726243699183035927701680494030_RKB, +0.191592375610517658429468812562799404E+0_RKB &
                                                                , +1.204166830144381512972585409228431810_RKB, +0.193215832729717172651871738633649002E+0_RKB &
                                                                , +1.197157747879441555149951169947666830_RKB, +0.194848709578749458096518050299059063E+0_RKB &
                                                                , +1.190125515426692069320479795685490130_RKB, +0.196491146544803628100426175257282935E+0_RKB &
                                                                , +1.183069142682686761029921497512774320_RKB, +0.198143288368357757303338740744073695E+0_RKB &
                                                                , +1.175987612015452098437940687032687740_RKB, +0.199805284331549509895377775340098004E+0_RKB &
                                                                , +1.168879876730833138384076802216189500_RKB, +0.201477288457465010601480272621186237E+0_RKB &
                                                                , +1.161744859445611442407678658730303350_RKB, +0.203159459721132113428367115004471988E+0_RKB &
                                                                , +1.154581450359927740740669252500180290_RKB, +0.204851962273072525723359017960643915E+0_RKB &
                                                                , +1.147388505420849058458757243582079270_RKB, +0.206554965676342511391397805599266881E+0_RKB &
                                                                , +1.140164844368151242664448937938463130_RKB, +0.208268645158074945659692042041697417E+0_RKB &
                                                                , +1.132909248652533753115898715892141530_RKB, +0.209993181876627252316799052585710853E+0_RKB &
                                                                , +1.125620459215533391227475593856507500_RKB, +0.211728763205541276296151308493989896E+0_RKB &
                                                                , +1.118297174119344981546693309349106980_RKB, +0.213475583035633628019879947303927958E+0_RKB &
                                                                , +1.110938046013575721417560779257849040_RKB, +0.215233842096659846410648766825208954E+0_RKB &
                                                                , +1.103541679424639718350709546845052450_RKB, +0.217003748300134423613077654918982009E+0_RKB &
                                                                , +1.096106627852021437094938365415991340_RKB, +0.218785517105043099072370007214575357E+0_RKB &
                                                                , +1.088631390653979821403993213445028160_RKB, +0.220579371908355906900178024260990122E+0_RKB &
                                                                , +1.081114409703403838077538136757234410_RKB, +0.222385544462441594762945169415386684E+0_RKB &
                                                                , +1.073554065792436288292686500464555760_RKB, +0.224204275321698924941276382479365999E+0_RKB &
                                                                , +1.065948674762122501754284026538321330_RKB, +0.226035814320961132564753113117134625E+0_RKB &
                                                                , +1.058296483330675084513966896229645790_RKB, +0.227880421088500051269826182271994471E+0_RKB &
                                                                , +1.050595664590929901514279329073961990_RKB, +0.229738365596760291802256320937473869E+0_RKB &
                                                                , +1.042844313144148970939980229020163640_RKB, +0.231609928754296215469589885980076311E+0_RKB &
                                                                , +1.035040439833440875887862746039381010_RKB, +0.233495403042770916853974440099008843E+0_RKB &
                                                                , +1.027181966035645772350639887939364020_RKB, +0.235395093203313594481836509609421923E+0_RKB &
                                                                , +1.019266717465484244962537467525869270_RKB, +0.237309316977027237519675921574483765E+0_RKB &
                                                                , +1.011292417439995739530489041191619980_RKB, +0.239238405905001516168070540916910782E+0_RKB &
                                                                , +1.003256679544672977063838476902952280_RKB, +0.241182706193826750302150947619295807E+0_RKB &
                                                                , +0.995156999635090923837912834290083799_RKB, +0.243142579653336370820347935962210089E+0_RKB &
                                                                , +0.986990747099062472368807733566367723_RKB, +0.245118404714142202924320107980890993E+0_RKB &
                                                                , +0.978755155294224603880824209605885949_RKB, +0.247110577533486783260721039414118216E+0_RKB &
                                                                , +0.970447311064224450680967868251783851_RKB, +0.249119513199040723215271175059742679E+0_RKB &
                                                                , +0.962064143223040583869351754884239000_RKB, +0.251145647041545878684196026994747375E+0_RKB &
                                                                , +0.953602409881086036147944393710646946_RKB, +0.253189436068676791163884011564958290E+0_RKB &
                                                                , +0.945058684468165463037907528952905452_RKB, +0.255251360534199633721355034811644515E+0_RKB &
                                                                , +0.936429340286575141234349646203595494_RKB, +0.257331925658493337610224432667231202E+0_RKB &
                                                                , +0.927710533402000123870677193352685174_RKB, +0.259431663518814566212074318832185428E+0_RKB &
                                                                , +0.918898183649590612180034442455129809_RKB, +0.261551135130401113878336401593219032E+0_RKB &
                                                                , +0.909987953496718494483529567690366580_RKB, +0.263690932742695838070411741690697980E+0_RKB &
                                                                , +0.900975224461221833746547856386022859_RKB, +0.265851682378732029260877257163356597E+0_RKB &
                                                                , +0.891855070732941566850586789359498648_RKB, +0.268034046650170417793435080776528378E+0_RKB &
                                                                , +0.882622229585165554772936671621819385_RKB, +0.270238727885765625479222906930855857E+0_RKB &
                                                                , +0.873271068088860754125762716220304264_RKB, +0.272466471617349972495381721975295691E+0_RKB &
                                                                , +0.863795545553308854813178505394360638_RKB, +0.274718070474985863152207934163143218E+0_RKB &
                                                                , +0.854189171008163807454180162989399104_RKB, +0.276994368552045124993231875042116081E+0_RKB &
                                                                , +0.844444954909153918889582440732428571_RKB, +0.279296266311992916126496729438627603E+0_RKB &
                                                                , +0.834555354086382178924726895193815695_RKB, +0.281624726122054644081592741705012005E+0_RKB &
                                                                , +0.824512208752292130518310689451625822_RKB, +0.283980778515329245279710753274268026E+0_RKB &
                                                                , +0.814306670135215230392694899997583155_RKB, +0.286365529303059635938328101021065686E+0_RKB &
                                                                , +0.803929116989971220407539518038181495_RKB, +0.288780167683694362995818015176624321E+0_RKB &
                                                                , +0.793369058840623296211094246674767659_RKB, +0.291225975526402318340687192537257840E+0_RKB &
                                                                , +0.782615023307233120893558043746138323_RKB, +0.293704338045591948359518681668325255E+0_RKB &
                                                                , +0.771654424224568084749572873233895398_RKB, +0.296216756132081268041022602880704688E+0_RKB &
                                                                , +0.760473406430108029348112105521972518_RKB, +0.298764860669019852359816171738518807E+0_RKB &
                                                                , +0.749056662017815292302576863549466866_RKB, +0.301350429240767224099410447613497339E+0_RKB &
                                                                , +0.737387211434295591278302895320305406_RKB, +0.303975405746574722328106158035369533E+0_RKB &
                                                                , +0.725446140909999639160421404681071701_RKB, +0.306641923566284096197801960631070354E+0_RKB &
                                                                , +0.713212285190975958395437234225351400_RKB, +0.309352333103853461802712328689370145E+0_RKB &
                                                                , +0.700661841106815072627797458756938969_RKB, +0.312109234772742761743862339129657007E+0_RKB &
                                                                , +0.687767892795788534294858623414951669_RKB, +0.314915518808718397997156048518493503E+0_RKB &
                                                                , +0.674499822837293822822291441953147501_RKB, +0.317774413735206587646816278462680542E+0_RKB &
                                                                , +0.660822574244419738417074112845703024_RKB, +0.320689545915737438541580666581786760E+0_RKB &
                                                                , +0.646695714894993817513389454402005337_RKB, +0.323665013485929872609853214812546138E+0_RKB &
                                                                , +0.632072236386061170945000136838048409_RKB, +0.326705479185682685503639905197287285E+0_RKB &
                                                                , +0.616896990007751449983468424580333684_RKB, +0.329816288403562001297942710464242773E+0_RKB &
                                                                , +0.601104617755992621533900682881952269_RKB, +0.333003621412417574590814956444079152E+0_RKB &
                                                                , +0.584616766106379321441587601292714165_RKB, +0.336274692838645610198728629890296241E+0_RKB &
                                                                , +0.567338257053818748196811566000406618_RKB, +0.339638017760732372141674530066184563E+0_RKB &
                                                                , +0.549151702327165120668100504067842453_RKB, +0.343103774061966311836681659618112095E+0_RKB &
                                                                , +0.529909720661558116786810165407173123_RKB, +0.346684307694080825876172229626686798E+0_RKB &
                                                                , +0.509423329602091814469823299066412422_RKB, +0.350394856987006815154163453112784578E+0_RKB &
                                                                , +0.487443966139236039301073245095571196_RKB, +0.354254625523424624005573055746854654E+0_RKB &
                                                                , +0.463634336790882217507922976793371502_RKB, +0.358288435101351989902170775658332694E+0_RKB &
                                                                , +0.437518402207871681933515025173607280_RKB, +0.362529398255472763687561092414814622E+0_RKB &
                                                                , +0.408389134611991145290558016416221589_RKB, +0.367023508970343116977068094647372845E+0_RKB &
                                                                , +0.375121332878380591495093443929722987_RKB, +0.371838172174307847618293715165674310E+0_RKB &
                                                                , +0.335737519214425235638195399148542308_RKB, +0.377079825919318504715611620323180171E+0_RKB &
                                                                , +0.286174591792072510002201653739714122_RKB, +0.382936353792390580381634933484297406E+0_RKB &
                                                                , +0.215241895984881699325976137068393403_RKB, +0.389807180887844207493515693365586429E+0_RKB &
                                                                , +0.000000000000000000000000000000000000_RKB, +0.398942280401432677939946059934381874E+0_RKB ], shape = [2, 128])

    !>  \brief
    !>  The constant array of type `real` of kind \RKB of shape `(1 : 2, 0 : 256)` containing the default 256-layers Ziggurat set
    !>  information that is used within [setNormRand](@ref pm_distNorm::setNormRand) to generate standard Normal random numbers.<br>
    !>
    !>  \details
    !>  The subset `ZIG_RKB(1, :)` corresponding to the lower right corners of the 256 Ziggurat rectangles are computed via
    !>  [getZigNorm](@ref pm_distNorm::getZigNorm) yielding a maximum absolute error `abserr = +0.562732655625645475814793555099033451E-33`
    !>  in the rectangle areas.<br>
    !>
    !>  \warning
    !>  This constant vector was generated on an `amd64` processor with the highest machine precision available for `real` type,
    !>  yielding a maximum of [ZIG_PRECISION](@ref pm_distNorm::ZIG_PRECISION) digits of precision.<br>
    !>  As such, care must be taken to not use this default Ziggurat set on processors with higher precision than used to generate this set.<br>
    !>
    !>  \see
    !>  [setZig](@ref pm_ziggurat::getZig)
    !>  [ZIG_RKB](@ref pm_distNorm::ZIG_RKB)
    !>  [setNormRand](@ref pm_distNorm::setNormRand)
    !>  [ZIG_PRECISION](@ref pm_distNorm::ZIG_PRECISION)
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{ZIG_RKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB), parameter :: ZIG_RKB(2, 0:256) = reshape([ZIGSET1, ZIGSET2], shape = [2, 257])

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar or array of arbitrary rank of random values from the univariate Normal distribution
    !>  with the specified input `mean` and optionally, with the specified input standard deviation `std` of the Normal distribution.
    !>
    !>  \details
    !>  The procedures of this generic interface are merely convenient wrappers
    !>  around the procedures of [setNormRand](@ref pm_distNorm::setNormRand).<br>
    !>
    !>  \param[in]      mean    :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `rand`, representing the mean of the Normal distribution.<br>
    !>  \param[in]      std     :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `rand`, representing the standard deviation of the Normal distribution.<br>
    !>                              (**optional**, default = `1.`, can be present only if `mean` is also present.)
    !>
    !>  \return
    !>  `rand`                  :   The output scalar or array of arbitrary rank, of<br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the Normal-distributed random output value.<br>
    !>
    !>  \interface{getNormRand}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormRand
    !>
    !>      rand = getNormRand(mean, std = std)
    !>      rand = getNormRand(rng, mean, std = std)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < std` must hold for the corresponding input arguments.<br>
    !>  All warnings and notes for [setNormRand](@ref pm_distNorm::setNormRand) also apply to the procedures of this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  On input, the procedures of this generic interface minimally require
    !>  the argument `mean` to allow compile-time resolution of the procedure calls.<br>
    !>
    !>  \remark
    !>  The procedures of this generic interface internally use the intrinsic Fortran RNG for uniform random number generation.<br>
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [setNormRandBox](@ref pm_distNorm::setNormRandBox)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [getNormCDF](@ref pm_distNorm::getNormCDF)<br>
    !>
    !>  \example{getNormRand}
    !>  \include{lineno} example/pm_distNorm/getNormRand/main.F90
    !>  \compilef{getNormRand}
    !>  \output{getNormRand}
    !>  \include{lineno} example/pm_distNorm/getNormRand/main.out.F90
    !>  \postproc{getNormRand}
    !>  \include{lineno} example/pm_distNorm/getNormRand/main.py
    !>  \vis{getNormRand}
    !>  \image html pm_distNorm/getNormRand/getNormRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)<br>
    !>
    !>  \final{getNormRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getNormRandRDMASA_D0_RK5(mean, std) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormRandRDMASA_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: mean
        real(RKG)   , intent(in)    , optional      :: std
        real(RKG)                                   :: rand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getNormRandRDMASA_D0_RK4(mean, std) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormRandRDMASA_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: mean
        real(RKG)   , intent(in)    , optional      :: std
        real(RKG)                                   :: rand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getNormRandRDMASA_D0_RK3(mean, std) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormRandRDMASA_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: mean
        real(RKG)   , intent(in)    , optional      :: std
        real(RKG)                                   :: rand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getNormRandRDMASA_D0_RK2(mean, std) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormRandRDMASA_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: mean
        real(RKG)   , intent(in)    , optional      :: std
        real(RKG)                                   :: rand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getNormRandRDMASA_D0_RK1(mean, std) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormRandRDMASA_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: mean
        real(RKG)   , intent(in)    , optional      :: std
        real(RKG)                                   :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar or array of arbitrary rank of random values from the standard univariate Normal distribution.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface use the Ziggurat rejection method to generate Normal-distributed random numbers.<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG for Gamma RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG for Normal RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type), implying the use of the intrinsic Fortran RNG.)
    !>  \param[out]     rand    :   The output
    !>                              <ol>
    !>                                  <li>    scalar, or<br>
    !>                                  <li>    array of rank `1`, or<br>
    !>                                  <li>    array of arbitrary rank if the `rng` argument is missing or set to [rngf_type](@ref pm_distUnif::rngf_type), or<br>
    !>                              </ol>
    !>                              of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL.
    !>                              </ol>
    !>                              On output, it contains standard Normal-distributed random value(s).<br>
    !>  \param[in]      zig     :   The input matrix of shape `(1:2, 0:*)` of the same type and kind as the output argument `rand`,
    !>                              containing the information about the Ziggurat layers.<br>
    !>                              This matrix is directly (and must have been) returned by [getZigNorm](@ref pm_distNorm::getZigNorm).<br>
    !>                              (**optional**. default = [ZIG_RKB](@ref pm_distNorm::ZIG_RKB))
    !>
    !>  \interface{setNormRand}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: setNormRand
    !>
    !>      call setNormRand(rand(..))
    !>      call setNormRand(rand(..), zig(:,:))
    !>
    !>      call setNormRand(rand(:))
    !>      call setNormRand(rand(:), zig(:,:))
    !>
    !>      call setNormRand(rng, rand)
    !>      call setNormRand(rng, rand, zig(:,:))
    !>
    !>      call setNormRand(rng, rand(:))
    !>      call setNormRand(rng, rand(:), zig(:,:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `precision(rand) <= ZIG_PRECISION` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0._RKG <= zig)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(zig, 1) == 2` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the input argument `rng` is missing.<br>
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `zig` is present
    !>  or `rng` is set to an object of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [setNormRandBox](@ref pm_distNorm::setNormRandBox)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [getNormCDF](@ref pm_distNorm::getNormCDF)<br>
    !>
    !>  \example{setNormRand}
    !>  \include{lineno} example/pm_distNorm/setNormRand/main.F90
    !>  \compilef{setNormRand}
    !>  \output{setNormRand}
    !>  \include{lineno} example/pm_distNorm/setNormRand/main.out.F90
    !>  \postproc{setNormRand}
    !>  \include{lineno} example/pm_distNorm/setNormRand/main.py
    !>  \vis{setNormRand}
    !>  \image html pm_distNorm/setNormRand/setNormRand.RK.png width=700
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{ziggurat_resolution, The runtime performance of [setNormRand](@ref pm_distNorm::setNormRand) for different implementations of the Box-Muller algorithm.}
    !>  \include{lineno} benchmark/pm_distNorm/ziggurat_resolution/main.F90
    !>  \compilefb{ziggurat_resolution}
    !>  \postprocb{ziggurat_resolution}
    !>  \include{lineno} benchmark/pm_distNorm/ziggurat_resolution/main.py
    !>  \visb{ziggurat_resolution}
    !>  \image html benchmark/pm_distNorm/ziggurat_resolution/benchmark.ziggurat_resolution.runtime.png width=1000
    !>  \image html benchmark/pm_distNorm/ziggurat_resolution/benchmark.ziggurat_resolution.runtime.ratio.png width=1000
    !>  \moralb{ziggurat_resolution}
    !>      -#  The benchmark procedures named `setNormRandX256` and `setNormRandFRNG` call the generic interface [setNormRand](@ref pm_distNorm::setNormRand)
    !>          to generate Normal random values via Ziggurat algorithm with the [Xoshiro256**](@ref pm_distUnif::xoshiro256ssw_type)
    !>          and the intrinsic Fortran Uniform RNG, respectively.<br>
    !>          This benchmark confirms the popular choice of the number of `256` Ziggurat layers as the optimal value.<br>
    !>          Nevertheless, the benchmarks appear to show improvements with even larger number of Ziggurat layers, depending on the platform, processor, and compiler choice.<br>
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{setNormRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setNormRandUDZD_D0_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setNormRandUDZD_D0_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setNormRandUDZD_D0_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setNormRandUDZD_D0_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setNormRandUDZD_D0_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setNormRandUFZD_D0_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setNormRandUFZD_D0_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setNormRandUFZD_D0_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setNormRandUFZD_D0_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setNormRandUFZD_D0_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setNormRandUXZD_D0_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setNormRandUXZD_D0_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setNormRandUXZD_D0_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setNormRandUXZD_D0_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setNormRandUXZD_D0_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setNormRandUDZA_D0_RK5(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setNormRandUDZA_D0_RK4(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setNormRandUDZA_D0_RK3(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setNormRandUDZA_D0_RK2(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setNormRandUDZA_D0_RK1(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setNormRandUFZA_D0_RK5(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setNormRandUFZA_D0_RK4(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setNormRandUFZA_D0_RK3(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setNormRandUFZA_D0_RK2(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setNormRandUFZA_D0_RK1(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setNormRandUXZA_D0_RK5(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setNormRandUXZA_D0_RK4(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setNormRandUXZA_D0_RK3(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setNormRandUXZA_D0_RK2(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setNormRandUXZA_D0_RK1(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setNormRandUDZD_D1_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setNormRandUDZD_D1_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setNormRandUDZD_D1_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setNormRandUDZD_D1_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setNormRandUDZD_D1_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setNormRandUFZD_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setNormRandUFZD_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setNormRandUFZD_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setNormRandUFZD_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setNormRandUFZD_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setNormRandUXZD_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setNormRandUXZD_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setNormRandUXZD_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setNormRandUXZD_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setNormRandUXZD_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setNormRandUDZA_D1_RK5(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setNormRandUDZA_D1_RK4(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setNormRandUDZA_D1_RK3(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setNormRandUDZA_D1_RK2(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setNormRandUDZA_D1_RK1(rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUDZA_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setNormRandUFZA_D1_RK5(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setNormRandUFZA_D1_RK4(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setNormRandUFZA_D1_RK3(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setNormRandUFZA_D1_RK2(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setNormRandUFZA_D1_RK1(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUFZA_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)             , intent(in)                    :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setNormRandUXZA_D1_RK5(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setNormRandUXZA_D1_RK4(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setNormRandUXZA_D1_RK3(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setNormRandUXZA_D1_RK2(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setNormRandUXZA_D1_RK1(rng, rand, zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandUXZA_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKG)                   , intent(out)                   :: rand(:)
        real(RKG)                   , intent(in)    , contiguous    :: zig(:, 0 :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar or array of arbitrary rank of random values from the univariate Normal distribution,
    !>  using the Box-Muller algorithm.
    !>
    !>  \details
    !>  The procedures of this generic interface use the Box-Muller transform method to convert
    !>  two uniformly distributed random input numbers to two output Normal-distributed random numbers.<br>
    !>  The procedures deliberately return two output random values to ensure the purity of the procedures.<br>
    !>  A more flexible implementation is available under the generic interface [getNormRand](@ref pm_distNorm::getNormRand).<br>
    !>
    !>  \param[inout]   rand1   :   The input/output scalar or array of the same shape as other array-like arguments, of,<br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL.
    !>                              </ol>
    !>                              On input, it must contain a uniformly-distributed random number in the range \f$[0, 1)\f$.<br>
    !>                              On output, it contains a Normal-distributed random value.<br>
    !>  \param[inout]   rand2   :   The input/output scalar or array of the same shape as other array-like arguments, of the same type and kind as `rand1`.<br>
    !>                              On input, it must contain a uniformly-distributed random number in the range \f$[0, 1)\f$.<br>
    !>                              On output, it contains a Normal-distributed random value.<br>
    !>  \param[out]     failed  :   The output scalar or array of the same shape as other array-like arguments, of type `logical` of default kind \LK
    !>                              that is `.true.` <b>if and only if</b> the Box-Muller rejection method fails to generate two Normal-random values.<br>
    !>                              In such a case, the procedure should be called until `failed = .false.` is returned on output indicating success.<br>
    !>                              (**optional**. If missing, the trigonometric Normal random number generation method will be used
    !>                              which is slightly (by about \f$\ms{25%}\f$) slower but guaranteed to succeed.)
    !>
    !>  \interface{setNormRandBox}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_distNorm, only: setNormRandBox
    !>      logical(LK) :: failed
    !>
    !>      ! Box-Muller: Trigonometric-method implementation.
    !>
    !>      call random_number(rand1)
    !>      call random_number(rand2)
    !>      call setNormRandBox(rand1, rand2)
    !>
    !>      ! Box-Muller: Rejection-method implementation.
    !>
    !>      do
    !>          call random_number(rand1)
    !>          call random_number(rand2)
    !>          call setNormRandBox(rand1, rand2, failed)
    !>          if (.not. failed) exit
    !>      end do
    !>
    !>  \endcode
    !>
    !>  The condition `0. <= rand1 .and. rand1 < 1.` must hold for the corresponding input arguments.<br>
    !>  The condition `0. <= rand2 .and. rand2 < 1.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [setNormRandBox](@ref pm_distNorm::setNormRandBox)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [getNormCDF](@ref pm_distNorm::getNormCDF)<br>
    !>  [setNormRandBox](@ref pm_distNorm::setNormRandBox)<br>
    !>
    !>  \example{setNormRandBox}
    !>  \include{lineno} example/pm_distNorm/setNormRandBox/main.F90
    !>  \compilef{setNormRandBox}
    !>  \output{setNormRandBox}
    !>  \include{lineno} example/pm_distNorm/setNormRandBox/main.out.F90
    !>  \postproc{setNormRandBox}
    !>  \include{lineno} example/pm_distNorm/setNormRandBox/main.py
    !>  \vis{setNormRandBox}
    !>  \image html pm_distNorm/setNormRandBox/setNormRandBox.RK.png width=700
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setNormRandBox_Basic_vs_Polar, The runtime performance of [setNormRandBox](@ref pm_distNorm::setNormRandBox) for different implementations of the Box-Muller algorithm.}
    !>  \include{lineno} benchmark/pm_distNorm/setNormRandBox_Basic_vs_Polar/main.F90
    !>  \compilefb{setNormRandBox_Basic_vs_Polar}
    !>  \postprocb{setNormRandBox_Basic_vs_Polar}
    !>  \include{lineno} benchmark/pm_distNorm/setNormRandBox_Basic_vs_Polar/main.py
    !>  \visb{setNormRandBox_Basic_vs_Polar}
    !>  \image html benchmark/pm_distNorm/setNormRandBox_Basic_vs_Polar/benchmark.setNormRandBox_Basic_vs_Polar.runtime.png width=1000
    !>  \image html benchmark/pm_distNorm/setNormRandBox_Basic_vs_Polar/benchmark.setNormRandBox_Basic_vs_Polar.runtime.ratio.png width=1000
    !>  \moralb{setNormRandBox_Basic_vs_Polar}
    !>      -#  The benchmark procedures named `setNormRandBoxBasic` and `setNormRandBoxPolar` call the generic interface [setNormRandBox](@ref pm_distNorm::setNormRandBox)
    !>          with the Basic and Polar (rejection sampling) implementations of the Box-Muller method, respectively.<br>
    !>          While both implementations are pure subroutines, the third benchmark implementation is `setNormRandBoxPolarImpure`
    !>          is an `impure` implementation of the Polar method that internally that does not require passing uniform random numbers to the subroutine.<br>
    !>      -#  The benchmark results confirm the widespread community observation that the Polar implementation of the Box-Muller method (which uses rejection sampling)
    !>          is slightly (by about 25%%) faster than the Basic implementation which involves two rather costly internal trigonometric function evaluations.<br>
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{setNormRandBox}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNormRandBox

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormRandBoxBasicDD_D0_RK5(rand1, rand2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormRandBoxBasicDD_D0_RK4(rand1, rand2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormRandBoxBasicDD_D0_RK3(rand1, rand2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormRandBoxBasicDD_D0_RK2(rand1, rand2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormRandBoxBasicDD_D0_RK1(rand1, rand2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNormRandBoxPolarDD_D0_RK5(rand1, rand2, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarDD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        logical(LK) , intent(out)                   :: failed
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNormRandBoxPolarDD_D0_RK4(rand1, rand2, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarDD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        logical(LK) , intent(out)                   :: failed
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNormRandBoxPolarDD_D0_RK3(rand1, rand2, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarDD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        logical(LK) , intent(out)                   :: failed
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNormRandBoxPolarDD_D0_RK2(rand1, rand2, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarDD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        logical(LK) , intent(out)                   :: failed
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNormRandBoxPolarDD_D0_RK1(rand1, rand2, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarDD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        logical(LK) , intent(out)                   :: failed
        real(RKG)   , intent(inout)                 :: rand1, rand2
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!
!#if RK5_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMD_RK5(rand1, rand2, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMD_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMD_RK4(rand1, rand2, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMD_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMD_RK3(rand1, rand2, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMD_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMD_RK2(rand1, rand2, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMD_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMD_RK1(rand1, rand2, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMD_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMS_RK5(rand1, rand2, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMS_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMS_RK4(rand1, rand2, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMS_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMS_RK3(rand1, rand2, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMS_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMS_RK2(rand1, rand2, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMS_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module subroutine setNormRandBoxBasicMS_RK1(rand1, rand2, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxBasicMS_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMD_RK5(rand1, rand2, failed, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMD_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMD_RK4(rand1, rand2, failed, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMD_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMD_RK3(rand1, rand2, failed, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMD_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMD_RK2(rand1, rand2, failed, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMD_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMD_RK1(rand1, rand2, failed, mean)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMD_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMS_RK5(rand1, rand2, failed, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMS_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMS_RK4(rand1, rand2, failed, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMS_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMS_RK3(rand1, rand2, failed, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMS_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMS_RK2(rand1, rand2, failed, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMS_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module subroutine setNormRandBoxPolarMS_RK1(rand1, rand2, failed, mean, std)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setNormRandBoxPolarMS_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        logical(LK) , intent(out)                   :: failed
!        real(RKG)   , intent(inout)                 :: rand1, rand2
!        real(RKG)   , intent(in)                    :: mean, std
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the lower right edges of the rectangles of a Ziggurat partitioning of the Normal density function
    !>  (and the corresponding density function values) to be used for Normal random number generation using the Ziggurat algorithm.
    !>
    !>  \details
    !>  See the documentation of [pm_ziggurat](@ref pm_ziggurat) for information on the Ziggurat method.<br>
    !>  See the documentation of [getZig](@ref pm_ziggurat::getZig) for the meaning of the elements of the output array.<br>
    !>
    !>  \param[in]  nlay    :   See the documentation of the corresponding argument of [getZig](@ref pm_ziggurat::getZig).<br>
    !>  \param[in]  abserr  :   See the documentation of the corresponding argument of [getZig](@ref pm_ziggurat::getZig).<br>
    !>  \param[in]  abstol  :   See the documentation of the corresponding argument of [getZig](@ref pm_ziggurat::getZig).<br>
    !>
    !>  \return
    !>  `zig`               :   See the documentation of the corresponding output of [getZig](@ref pm_ziggurat::getZig).<br>
    !>
    !>  \interface{getZigNorm}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getZigNorm
    !>
    !>      zig(1 : 2, 0 : nlay) = getZigNorm(nlay, abserr, abstol = abstol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All conditions that hold for [getZig](@ref pm_ziggurat::getZig) equally apply to the procedures of this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getZig](@ref pm_ziggurat::getZig)<br>
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>
    !>  \example{getZigNorm}
    !>  \include{lineno} example/pm_distNorm/getZigNorm/main.F90
    !>  \compilef{getZigNorm}
    !>  \output{getZigNorm}
    !>  \include{lineno} example/pm_distNorm/getZigNorm/main.out.F90
    !>  \postproc{getZigNorm}
    !>  \include{lineno} example/pm_distNorm/getZigNorm/main.py
    !>  \vis{getZigNorm}
    !>  \image html pm_distNorm/getZigNorm/getZigNorm.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)<br>
    !>
    !>  \final{getZigNorm}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getZigNorm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getZigNorm_RK5(nlay, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZigNorm_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: nlay
        real(RKG)   , intent(out)                   :: abserr
        real(RKG)   , intent(in)    , optional      :: abstol
        real(RKG)                                   :: zig(2, 0 : nlay)
    end function
#endif

#if RK4_ENABLED
    impure module function getZigNorm_RK4(nlay, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZigNorm_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: nlay
        real(RKG)   , intent(out)                   :: abserr
        real(RKG)   , intent(in)    , optional      :: abstol
        real(RKG)                                   :: zig(2, 0 : nlay)
    end function
#endif

#if RK3_ENABLED
    impure module function getZigNorm_RK3(nlay, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZigNorm_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: nlay
        real(RKG)   , intent(out)                   :: abserr
        real(RKG)   , intent(in)    , optional      :: abstol
        real(RKG)                                   :: zig(2, 0 : nlay)
    end function
#endif

#if RK2_ENABLED
    impure module function getZigNorm_RK2(nlay, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZigNorm_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: nlay
        real(RKG)   , intent(out)                   :: abserr
        real(RKG)   , intent(in)    , optional      :: abstol
        real(RKG)                                   :: zig(2, 0 : nlay)
    end function
#endif

#if RK1_ENABLED
    impure module function getZigNorm_RK1(nlay, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZigNorm_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: nlay
        real(RKG)   , intent(out)                   :: abserr
        real(RKG)   , intent(in)    , optional      :: abstol
        real(RKG)                                   :: zig(2, 0 : nlay)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \cond excluded
!    !   Internal module function used for the computation of the ziggurat partitions.
!    interface getFuncNorm
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    pure elemental module function getFuncNorm_RK5(x) result(func)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncNorm_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: func
!    end function
!#endif
!
!#if RK4_ENABLED
!    pure elemental module function getFuncNorm_RK4(x) result(func)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncNorm_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: func
!    end function
!#endif
!
!#if RK3_ENABLED
!    pure elemental module function getFuncNorm_RK3(x) result(func)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncNorm_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: func
!    end function
!#endif
!
!#if RK2_ENABLED
!    pure elemental module function getFuncNorm_RK2(x) result(func)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncNorm_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: func
!    end function
!#endif
!
!#if RK1_ENABLED
!    pure elemental module function getFuncNorm_RK1(x) result(func)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncNorm_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: func
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!    !>  \endcond excluded
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \cond excluded
!    !   Internal module function used for the computation of the ziggurat partitions.
!    interface getGradNorm
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    pure elemental module function getGradNorm_RK5(x) result(grad)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGradNorm_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: grad
!    end function
!#endif
!
!#if RK4_ENABLED
!    pure elemental module function getGradNorm_RK4(x) result(grad)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGradNorm_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: grad
!    end function
!#endif
!
!#if RK3_ENABLED
!    pure elemental module function getGradNorm_RK3(x) result(grad)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGradNorm_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: grad
!    end function
!#endif
!
!#if RK2_ENABLED
!    pure elemental module function getGradNorm_RK2(x) result(grad)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGradNorm_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: grad
!    end function
!#endif
!
!#if RK1_ENABLED
!    pure elemental module function getGradNorm_RK1(x) result(grad)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGradNorm_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in)    :: x
!        real(RKG)                   :: grad
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!    !>  \endcond excluded
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \cond excluded
!    !   Internal module function used for the computation of the ziggurat partitions.
!    interface getFuncInvNorm
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    pure elemental module function getFuncInvNorm_RK5(func) result(x)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncInvNorm_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)    :: func
!        real(RKG)                   :: x
!    end function
!#endif
!
!#if RK4_ENABLED
!    pure elemental module function getFuncInvNorm_RK4(func) result(x)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncInvNorm_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)    :: func
!        real(RKG)                   :: x
!    end function
!#endif
!
!#if RK3_ENABLED
!    pure elemental module function getFuncInvNorm_RK3(func) result(x)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncInvNorm_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in)    :: func
!        real(RKG)                   :: x
!    end function
!#endif
!
!#if RK2_ENABLED
!    pure elemental module function getFuncInvNorm_RK2(func) result(x)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncInvNorm_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in)    :: func
!        real(RKG)                   :: x
!    end function
!#endif
!
!#if RK1_ENABLED
!    pure elemental module function getFuncInvNorm_RK1(func) result(x)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getFuncInvNorm_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in)    :: func
!        real(RKG)                   :: x
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!    !>  \endcond excluded
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \cond excluded
!    !   Internal module function used for the computation of the ziggurat partitions.
!    interface getZigAreaNorm
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getZigAreaNorm_RK5(r) result(area)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigAreaNorm_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)    :: r
!        real(RKG)                   :: area
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getZigAreaNorm_RK4(r) result(area)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigAreaNorm_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)    :: r
!        real(RKG)                   :: area
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getZigAreaNorm_RK3(r) result(area)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigAreaNorm_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in)    :: r
!        real(RKG)                   :: area
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getZigAreaNorm_RK2(r) result(area)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigAreaNorm_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in)    :: r
!        real(RKG)                   :: area
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getZigAreaNorm_RK1(r) result(area)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigAreaNorm_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in)    :: r
!        real(RKG)                   :: area
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!    !>  \endcond excluded
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the entropy of the Normal distribution with the input natural logarithm of the variance.
    !>
    !>  \details
    !>  The entropy of the Normal distribution is defined by the following equation,
    !>  \f{equation}{
    !>      \large
    !>      \mathcal{H}(\sigma^2) = \frac{1}{2} \log(2\pi\sigma^2) + \frac{1}{2}
    !>  \f}
    !>
    !>  \param[in]  logVar  :   The input scalar or array of arbitrary rank and shape,
    !>                          of the same type and kind as the output `entropy`, containing the
    !>                          natural logarithm of the variance of the Normal distribution.<br>
    !>
    !>  \return
    !>  `entropy`           :   The output scalar or array of the same rank and shape as the input `logVar`, of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the entropy of the Normal distribution.
    !>
    !>  \interface{getNormEntropy}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormEntropy
    !>
    !>      entropy = getNormEntropy(logVar)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \example{getNormEntropy}
    !>  \include{lineno} example/pm_distNorm/getNormEntropy/main.F90
    !>  \compilef{getNormEntropy}
    !>  \output{getNormEntropy}
    !>  \include{lineno} example/pm_distNorm/getNormEntropy/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{getNormEntropy}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormEntropy

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getNormEntropy_RK5(logVar) result(entropy)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormEntropy_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                :: logVar
        real(RKG)                                           :: entropy
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getNormEntropy_RK4(logVar) result(entropy)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormEntropy_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                :: logVar
        real(RKG)                                           :: entropy
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getNormEntropy_RK3(logVar) result(entropy)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormEntropy_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                :: logVar
        real(RKG)                                           :: entropy
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getNormEntropy_RK2(logVar) result(entropy)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormEntropy_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                :: logVar
        real(RKG)                                           :: entropy
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getNormEntropy_RK1(logVar) result(entropy)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormEntropy_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                :: logVar
        real(RKG)                                           :: entropy
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormEntropy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the Fisher Information of the Normal distribution.
    !>
    !>  \details
    !>  The Fisher information for the Normal distribution is defined by the following equation,
    !>  \f{equation}{
    !>      \large
    !>      \mathcal{I}(\mu,\sigma) =
    !>      \begin{pmatrix}
    !>          \frac{1}{\sigma^2} & 0 \\
    !>          0 & \frac{2}{\sigma^2} \\
    !>      \end{pmatrix}
    !>  \f}
    !>
    !>  \param[in]  varInv  :   The input non-negative scalar of the same type and kind as the output `Fisher`
    !>                          representing the inverse of the variance of the Normal distribution.<br>
    !>
    !>  \return
    !>  `Fisher`            :   The output scalar of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL ,
    !>                          </ol>
    !>                          containing the Fisher Information matrix of the Normal distribution.
    !>
    !>  \interface{getNormFisher}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormFisher
    !>
    !>      Fisher(1:2,1:2) = getNormFisher(varInv)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. <= varInv` must hold for the correspondddding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \example{getNormFisher}
    !>  \include{lineno} example/pm_distNorm/getNormFisher/main.F90
    !>  \compilef{getNormFisher}
    !>  \output{getNormFisher}
    !>  \include{lineno} example/pm_distNorm/getNormFisher/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \final{getNormFisher}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormFisher

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormFisher_RK5(varInv) result(Fisher)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormFisher_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                :: varInv
        real(RKG)                                           :: Fisher(2,2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormFisher_RK4(varInv) result(Fisher)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormFisher_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                :: varInv
        real(RKG)                                           :: Fisher(2,2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormFisher_RK3(varInv) result(Fisher)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormFisher_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                :: varInv
        real(RKG)                                           :: Fisher(2,2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormFisher_RK2(varInv) result(Fisher)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormFisher_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                :: varInv
        real(RKG)                                           :: Fisher(2,2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormFisher_RK1(varInv) result(Fisher)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormFisher_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                :: varInv
        real(RKG)                                           :: Fisher(2,2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormFisher

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Kullback-Leibler Divergence (**KLD**) \f$D_{KL}(P \parallel Q)\f$
    !>  of a given univariate Normal distribution \f$Q\f$ from a reference Normal distribution \f$P\f$.
    !>
    !>  \details
    !>  The Kullback-Leibler Divergence, also known as the **relative entropy**, of a univariate Normal
    !>  distribution \f$Q\f$ from a reference univariate Normal distribution \f$P\f$ is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_{KL}(P \parallel Q) =
    !>      \frac{(\mu_P - \mu_Q)^2}{2\sigma_Q^2} +
    !>      \frac{1}{2}\bigg( \frac{\sigma_P}{\sigma_Q} - \ln\big(\frac{\sigma_P}{\sigma_Q}\big) - 1 \bigg) ~,
    !>  \f}
    !>
    !>  where \f$\mu\f$ and \f$\sigma\f$ represent the respective standard deviations of the Normal distributions.
    !>
    !>  \param[in]  meanDiffSq  :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              kind as the output `kld`, representing the difference-squared of the
    !>                              mean parameters of the two Normal distributions (\f$(P, Q)\f$).<br>
    !>                              (**optional**, default = `0.`, that is, the two distributions have the same mean.)
    !>  \param[in]  varP        :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as the output `kld` representing the variance of the **reference** Normal distribution \f$P\f$.<br>
    !>                              (**optional**, default = `1.`, it must be present <b>if and only if</b> `varQ` is also present.)
    !>  \param[in]  varQ        :   The input scalar or array of the same shape as other array-like arguments, of the same type and
    !>                              kind as `meanDiffSq` representing the variance of the **target** Normal distribution \f$Q\f$.<br>
    !>                              (**optional**, default = `1.`, it must be present <b>if and only if</b> `varP` is also present.)
    !>
    !>  \return
    !>  `kld`                   :   The input scalar or array of the same shape as input array-like arguments, of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              representing the Kullback-Leibler divergence of the target distribution with respect to the reference.
    !>
    !>  \interface{getNormKLD}
    !>  \code{.F90}
    !>
    !>      use pm_distNorm, only: getNormKLD
    !>
    !>      KLD = getNormKLD(meanDiffSq)
    !>      KLD = getNormKLD(varP, varQ)
    !>      KLD = getNormKLD(meanDiffSq, varP, varQ)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 <= meanDiffSq` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < varP` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 < varQ` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \example{getNormKLD}
    !>  \include{lineno} example/pm_distNorm/getNormKLD/main.F90
    !>  \compilef{getNormKLD}
    !>  \output{getNormKLD}
    !>  \include{lineno} example/pm_distNorm/getNormKLD/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distNorm](@ref test_pm_distNorm)
    !>
    !>  \todo
    !>  \phigh The KLD implementation should be implemented with `log1mexp` for better accuracy.
    !>
    !>  \final{getNormKLD}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormKLD

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNormKLDMD_RK5(meanDiffSq) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)                                           :: kld
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNormKLDMD_RK4(meanDiffSq) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)                                           :: kld
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNormKLDMD_RK3(meanDiffSq) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)                                           :: kld
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNormKLDMD_RK2(meanDiffSq) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)                                           :: kld
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNormKLDMD_RK1(meanDiffSq) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)                                           :: kld
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNormKLDDV_RK5(varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDDV_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNormKLDDV_RK4(varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDDV_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNormKLDDV_RK3(varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDDV_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNormKLDDV_RK2(varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDDV_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNormKLDDV_RK1(varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDDV_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNormKLDMV_RK5(meanDiffSq, varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMV_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNormKLDMV_RK4(meanDiffSq, varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMV_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNormKLDMV_RK3(meanDiffSq, varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMV_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNormKLDMV_RK2(meanDiffSq, varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMV_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNormKLDMV_RK1(meanDiffSq, varP, varQ) result(kld)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormKLDMV_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)                :: meanDiffSq
        real(RKG)               , intent(in)                :: varP, varQ
        real(RKG)                                           :: kld
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormKLD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distNorm
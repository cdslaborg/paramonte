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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Poisson distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Poisson distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Mass Function (**PMF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The Poisson distribution is a discrete probability distribution that
!>  expresses the probability of a given number of events occurring in a fixed interval of time or space if these events
!>  occur with a known constant mean rate and independently of the time since the last event.<br>
!>  It is named after French mathematician Siméon Denis Poisson.<br>
!>  The Poisson distribution can also be used for the number of events in other specified interval types such as distance, area, or volume.<br>
!>  A discrete random variable \f$X\f$ is said to have a Poisson distribution, with parameter \f$\lambda > 0\f$
!>  if it has a probability mass function given by,
!>  \f{equation}{
!>      f(k; \lambda) = \pi(X = k) = \frac {\lambda^k \exp\left(-\lambda\right)}{k!} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$k\f$ is the number of occurrences (\f$k = 0, 1, 2, \ldots\f$), and
!>      <li>    \f$\ms{!}\f$ is the factorial function.
!>  </ol>
!>  The positive real number \f$\lambda\f$ is equal to the expected value of \f$X\f$ and also to its variance.<br>
!>  \f{equation}{
!>      \lambda = \up{E}(X) = \up{Var}(X) ~.
!>  \f}
!>
!>  The **CDF** of the **Poisson distribution** with parameter \f$\lambda\f$ is defined as,
!>  \f{eqnarray}{
!>  \ms{CDF}(k | \lambda)
!>  &=& \frac{\Gamma(\lfloor k + 1 \rfloor, \lambda)}{\lfloor k \rfloor !} ~, \nonumber \\
!>  &=& \exp\left(-\lambda\right) \sum _{j=0}^{\lfloor k \rfloor}{\frac{\lambda^{j}}{j!}} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$k\f$ is the number of occurrences (\f$k = 0, 1, 2, \ldots\f$), and
!>      <li>    \f$\ms{!}\f$ is the factorial function, and
!>      <li>    \f$\Gamma(x, y) / \lfloor k \rfloor !\f$ is the [regularized upper incomplete gamma function](@ref pm_mathGamma::getGammaIncUpp), and
!>      <li>    \f$\lfloor k \rfloor\f$ is the [floor function](https://en.wikipedia.org/wiki/Floor_and_ceiling_functions).
!>  </ol>
!>
!>  **Random Number Generation**<br>
!>
!>  The RNG generic interfaces of this module use two different approaches for Poisson RNG
!>  for different ranges of \f$\lambda\f$ parameter values of the Poisson distribution.<br>
!>  <ol>
!>      <li>    When \f$\lambda <\f$ [LAMBDA_LIMIT](@ref pm_distPois::LAMBDA_LIMIT), a RNG algorithm due to
!>              [Donald Ervin Knuth](https://en.wikipedia.org/wiki/Donald_Knuth) is used.<br>
!>      <li>    When [LAMBDA_LIMIT](@ref pm_distPois::LAMBDA_LIMIT) \f$\leq \lambda\f$, a rejection-based RNG algorithm
!>              due to *Hormann, 1993, The transformed rejection method for generating Poisson random variables* is used.<br>
!>              This rejection method has an efficiency slight less than \f$90\%\f$.<br>
!>  </ol>
!>
!>  \see
!>  [pm_distPois](@ref pm_distPois)<br>
!>  [pm_distPois](@ref pm_distPois)<br>
!>  [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getPoisLogPMF_vs_setPoisLogPMF, The runtime performance of [getPoisLogPMF](@ref pm_distPois::getPoisLogPMF) vs. [setPoisLogPMF](@ref pm_distPois::setPoisLogPMF)}
!>  \include{lineno} benchmark/pm_distPois/getPoisLogPMF_vs_setPoisLogPMF/main.F90
!>  \compilefb{getPoisLogPMF_vs_setPoisLogPMF}
!>  \postprocb{getPoisLogPMF_vs_setPoisLogPMF}
!>  \include{lineno} benchmark/pm_distPois/getPoisLogPMF_vs_setPoisLogPMF/main.py
!>  \visb{getPoisLogPMF_vs_setPoisLogPMF}
!>  \image html benchmark/pm_distPois/getPoisLogPMF_vs_setPoisLogPMF/benchmark.getPoisLogPMF_vs_setPoisLogPMF.runtime.png width=1000
!>  \image html benchmark/pm_distPois/getPoisLogPMF_vs_setPoisLogPMF/benchmark.getPoisLogPMF_vs_setPoisLogPMF.runtime.ratio.png width=1000
!>  \moralb{getPoisLogPMF_vs_setPoisLogPMF}
!>      -#  The procedures under the generic interface [getPoisLogPMF](@ref pm_distPois::getPoisLogPMF) are functions while
!>          the procedures under the generic interface [setPoisLogPMF](@ref pm_distPois::setPoisLogPMF) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than
!>          the subroutine interface when the input `array` size is small.<br>
!>          Otherwise, the difference appears to be marginal and insignificant in most practical situations.<br>
!>
!>  \benchmark{setPoisLogPMF-logLambda-missing_vs_present, The runtime performance of [setPoisLogPMF](@ref pm_distPois::setPoisLogPMF) with and without `logLambda`}
!>  \include{lineno} benchmark/pm_distPois/setPoisLogPMF-logLambda-missing_vs_present/main.F90
!>  \compilefb{setPoisLogPMF-logLambda-missing_vs_present}
!>  \postprocb{setPoisLogPMF-logLambda-missing_vs_present}
!>  \include{lineno} benchmark/pm_distPois/setPoisLogPMF-logLambda-missing_vs_present/main.py
!>  \visb{setPoisLogPMF-logLambda-missing_vs_present}
!>  \image html benchmark/pm_distPois/setPoisLogPMF-logLambda-missing_vs_present/benchmark.setPoisLogPMF-logLambda-missing_vs_present.runtime.png width=1000
!>  \image html benchmark/pm_distPois/setPoisLogPMF-logLambda-missing_vs_present/benchmark.setPoisLogPMF-logLambda-missing_vs_present.runtime.ratio.png width=1000
!>  \moralb{setPoisLogPMF-logLambda-missing_vs_present}
!>      -#  The procedures under the generic interface [setPoisLogPMF](@ref pm_distPois::setPoisLogPMF)
!>          accept an extra argument `logLambda = log(lambda)` while the procedures under the generic interface
!>          [getPoisLogPMF](@ref pm_distPois::getPoisLogPMF) compute this term internally with every procedure call.<br>
!>          In the presence of this argument, the logarithmic computation `log(lambda)` will be avoided.<br>
!>          As such, the presence of `logLambda` is expected to lead to faster computations.<br>
!>
!>  \test
!>  [test_pm_distPois](@ref test_pm_distPois)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distPois

    use pm_kind, only: SK, IK, LK, RKB
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distPois"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Poisson
    !>  as defined in the description of [pm_distPois](@ref pm_distPois).
    !>
    !>  \details
    !>  See the documentation of [pm_distPois](@ref pm_distPois) for the definition of the Poisson distribution.
    !>
    !>  \interface{distPois_type}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: distPois_type
    !>      type(distPois_type) :: distPois
    !>
    !>      distPois = distPois_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distPois_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distPois_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Poisson distribution</b> for an input `count` within the discrete integer support of the distribution \f$[0, +\infty)\f$.
    !>
    !>  \param[in]  count       :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values at which the PMF must be computed.<br>
    !>  \param[in]  lambda      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the location parameter of the distribution.<br>
    !>
    !>  \return
    !>  `logPMF`                :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `lambda`, containing the natural logarithm of the PMF of
    !>                              the distribution at the specified point.<br>
    !>
    !>  \interface{getPoisLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: getPoisLogPMF
    !>
    !>      logPMF = getPoisLogPMF(count, lambda)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= count` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lambda` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPoisLogPMF](@ref pm_distPois::setPoisLogPMF)<br>
    !>
    !>  \example{getPoisLogPMF}
    !>  \include{lineno} example/pm_distPois/getPoisLogPMF/main.F90
    !>  \compilef{getPoisLogPMF}
    !>  \output{getPoisLogPMF}
    !>  \include{lineno} example/pm_distPois/getPoisLogPMF/main.out.F90
    !>  \postproc{getPoisLogPMF}
    !>  \include{lineno} example/pm_distPois/getPoisLogPMF/main.py
    !>  \vis{getPoisLogPMF}
    !>  \image html pm_distPois/getPoisLogPMF/getPoisLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \final{getPoisLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPoisLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPoisLogPMF_RK5(count, lambda) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisLogPMF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPoisLogPMF_RK4(count, lambda) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisLogPMF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPoisLogPMF_RK3(count, lambda) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisLogPMF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPoisLogPMF_RK2(count, lambda) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisLogPMF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPoisLogPMF_RK1(count, lambda) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisLogPMF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: logPMF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Poisson distribution</b> for an input `count` within the discrete integer support of the distribution \f$[0, +\infty)\f$.
    !>
    !>  \param[out] logPMF      :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `lambda`,
    !>                              containing the natural logarithm of the PMF of the distribution at the specified `count`.<br>
    !>  \param[in]  count       :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values at which the PMF must be computed.<br>
    !>  \param[in]  lambda      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the location parameter of the distribution.<br>
    !>  \param[in]  logLambda   :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `lambda`
    !>                              containing the natural logarithm of the location parameter of the PMF `lambda`.<br>
    !>                              Specifying this argument leads to significantly faster computations of the PMF.<br>
    !>                              (**optional**, default = `log(lambda)`.)
    !>
    !>  \interface{setPoisLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: setPoisLogPMF
    !>
    !>      call setPoisLogPMF(logPMF, count, lambda)
    !>      call setPoisLogPMF(logPMF, count, lambda, logLambda)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= count` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lambda` must hold for the corresponding input arguments.<br>
    !>  The condition `logLambda = log(lambda)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPoisLogPMF](@ref pm_distPois::getPoisLogPMF)<br>
    !>
    !>  \example{setPoisLogPMF}
    !>  \include{lineno} example/pm_distPois/setPoisLogPMF/main.F90
    !>  \compilef{setPoisLogPMF}
    !>  \output{setPoisLogPMF}
    !>  \include{lineno} example/pm_distPois/setPoisLogPMF/main.out.F90
    !>  \postproc{setPoisLogPMF}
    !>  \include{lineno} example/pm_distPois/setPoisLogPMF/main.py
    !>  \vis{setPoisLogPMF}
    !>  \image html pm_distPois/setPoisLogPMF/setPoisLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setPoisLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPoisLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPoisLogPMFDef_RK5(logPMF, count, lambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPoisLogPMFDef_RK4(logPMF, count, lambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPoisLogPMFDef_RK3(logPMF, count, lambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPoisLogPMFDef_RK2(logPMF, count, lambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPoisLogPMFDef_RK1(logPMF, count, lambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPoisLogPMFLog_RK5(logPMF, count, lambda, logLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFLog_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: lambda, logLambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPoisLogPMFLog_RK4(logPMF, count, lambda, logLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFLog_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: lambda, logLambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPoisLogPMFLog_RK3(logPMF, count, lambda, logLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFLog_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: lambda, logLambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPoisLogPMFLog_RK2(logPMF, count, lambda, logLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFLog_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: lambda, logLambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPoisLogPMFLog_RK1(logPMF, count, lambda, logLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisLogPMFLog_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: lambda, logLambda
        real(RKG)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  <b>Poisson distribution</b> for an input `count` within the discrete integer support of the distribution \f$[0, +\infty)\f$.
    !>
    !>  \param[in]  count       :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values at which the CDF must be computed.<br>
    !>  \param[in]  lambda      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the location parameter of the distribution.<br>
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `lambda`, containing the CDF of
    !>                              the distribution at the specified point.<br>
    !>
    !>  \interface{getPoisCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: getPoisCDF
    !>
    !>      cdf = getPoisCDF(count, lambda)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= count` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lambda` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPoisCDF](@ref pm_distPois::setPoisCDF)<br>
    !>
    !>  \example{getPoisCDF}
    !>  \include{lineno} example/pm_distPois/getPoisCDF/main.F90
    !>  \compilef{getPoisCDF}
    !>  \output{getPoisCDF}
    !>  \include{lineno} example/pm_distPois/getPoisCDF/main.out.F90
    !>  \postproc{getPoisCDF}
    !>  \include{lineno} example/pm_distPois/getPoisCDF/main.py
    !>  \vis{getPoisCDF}
    !>  \image html pm_distPois/getPoisCDF/getPoisCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \final{getPoisCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPoisCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getPoisCDF_RK5(count, lambda) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getPoisCDF_RK4(count, lambda) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getPoisCDF_RK3(count, lambda) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getPoisCDF_RK2(count, lambda) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getPoisCDF_RK1(count, lambda) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: count
        real(RKG)   , intent(in)                    :: lambda
        real(RKG)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the <b>Poisson distribution</b>.
    !>
    !>  \param[out] cdf             :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `lambda`,
    !>                                  containing the CDF of the distribution at the specified input value(s).<br>
    !>  \param[in]  countP1         :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as the input argument `lambda`, containing the value(s) at which the CDF must be computed plus `1.`.<br>
    !>                                  Although `countP1` is of type `real`, it must be a whole number.<br>
    !>                                  The enforcement of `real` type is to ensure the performance of this generic interface.<br>
    !>  \param[in]  lambda          :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the location parameter of the distribution.<br>
    !>  \param[out] info            :   The output scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm.<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the CDF**.<br>
    !>                                  For more information, see the documentation of [setGammaInc](@ref pm_mathGamma::setGammaInc).<br>
    !>
    !>  \interface{setPoisCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: setPoisCDF
    !>
    !>      call setPoisCDF(cdf, countP1, lambda, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lambda` must hold for the corresponding input arguments.<br>
    !>  The condition `1. <= countP1` must hold for the corresponding input arguments.<br>
    !>  The condition `mod(countP1, 1._RKG) == 0._RKG` must hold for the corresponding input arguments,
    !>  that is, the input argument `countP1` must be a whole number.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPoisCDF](@ref pm_distPois::getPoisCDF)<br>
    !>
    !>  \example{setPoisCDF}
    !>  \include{lineno} example/pm_distPois/setPoisCDF/main.F90
    !>  \compilef{setPoisCDF}
    !>  \output{setPoisCDF}
    !>  \include{lineno} example/pm_distPois/setPoisCDF/main.out.F90
    !>  \postproc{setPoisCDF}
    !>  \include{lineno} example/pm_distPois/setPoisCDF/main.py
    !>  \vis{setPoisCDF}
    !>  \image html pm_distPois/setPoisCDF/setPoisCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \final{setPoisCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPoisCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPoisCDFLog_RK5(cdf, countP1, lambda, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisCDFLog_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: countP1, lambda
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPoisCDFLog_RK4(cdf, countP1, lambda, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisCDFLog_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: countP1, lambda
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPoisCDFLog_RK3(cdf, countP1, lambda, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisCDFLog_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: countP1, lambda
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPoisCDFLog_RK2(cdf, countP1, lambda, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisCDFLog_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: countP1, lambda
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPoisCDFLog_RK1(cdf, countP1, lambda, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisCDFLog_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: countP1, lambda
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of kind \RKB, representing the value of the parameter of the Poisson distribution
    !>  above which the rejection method of *Hormann, 1993, The transformed rejection method for generating Poisson random variables*
    !>  for generating Poisson-distributed random values is valid.<br>
    !>
    !>  \details
    !>  This constant exists merely as a reference with which decision can be made about
    !>  the proper [setPoisRand](@ref pm_distPois::setPoisRand) interface usage.<br>
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \final{getPoisRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)   , parameter :: LAMBDA_LIMIT = 10._RKB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank of) random value(s) from the Poisson distribution.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPois](@ref pm_distPois) for more details.<br>
    !>
    !>  \param[in]  lambda  :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                          <ul>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          representing the location/scale parameter of the distribution.<br>
    !>
    !>  \return
    !>  `rand`              :   The output positive scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                          of type `integer` of default kind \IK, containing random value(s) from the distribution.<br>
    !>
    !>  \interface{getPoisRand}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: getPoisRand
    !>      use pm_kind, only: IK
    !>      integer(IK) :: rand
    !>
    !>      rand = getPoisRand(lambda)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lambda` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setPoisRand](@ref pm_distPois::setPoisRand)<br>
    !>
    !>  \example{getPoisRand}
    !>  \include{lineno} example/pm_distPois/getPoisRand/main.F90
    !>  \compilef{getPoisRand}
    !>  \output{getPoisRand}
    !>  \include{lineno} example/pm_distPois/getPoisRand/main.out.F90
    !>  \postproc{getPoisRand}
    !>  \include{lineno} example/pm_distPois/getPoisRand/main.py
    !>  \vis{getPoisRand}
    !>  \image html pm_distPois/getPoisRand/getPoisRand.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \final{getPoisRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getPoisRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getPoisRand_RK5(lambda) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisRand_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                                 :: rand
        real(RKG)   , intent(in)                    :: lambda
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getPoisRand_RK4(lambda) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisRand_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                                 :: rand
        real(RKG)   , intent(in)                    :: lambda
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getPoisRand_RK3(lambda) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisRand_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                                 :: rand
        real(RKG)   , intent(in)                    :: lambda
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getPoisRand_RK2(lambda) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisRand_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                                 :: rand
        real(RKG)   , intent(in)                    :: lambda
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getPoisRand_RK1(lambda) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPoisRand_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                                 :: rand
        real(RKG)   , intent(in)                    :: lambda
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Poisson distribution.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distPois](@ref pm_distPois) for more details.<br>
    !>
    !>  \param[inout]   rng             :   The input/output scalar that can be an object of,
    !>                                      <ol>
    !>                                          <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                                  implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                          <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                                  implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                                      </ol>
    !>                                      (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type), implying the use of the intrinsic Fortran URNG.)
    !>  \param[out]     rand            :   The output scalar or
    !>                                      <ol>
    !>                                          <li>    array of rank `1`, or<br>
    !>                                          <li>    array of arbitrary rank if the `rng` argument is missing or set to [rngf_type](@ref pm_distUnif::rngf_type),<br>
    !>                                      </ol>
    !>                                      of type `integer` of default kind \IK, containing the Poisson-distributed random value(s).<br>
    !>  \param[in]      lambda          :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                      <ol>
    !>                                          <li>    type `real` of kind \RKALL,
    !>                                      </ol>
    !>                                      representing the parameter(s) of the distribution(s).<br>
    !>                                      (**optional**. It must be present **if and only if** [LAMBDA_LIMIT](@ref pm_distPois::LAMBDA_LIMIT) \f$\leq \lambda\f$)
    !>  \param[in]      logLambda       :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                      of the same type and kind as `lambda`, representing `log(lambda)`.<br>
    !>                                      This argument is required to ensure fast computation of many random values from the same distribution.<br>
    !>                                      (**optional**. It must be present **if and only if** [LAMBDA_LIMIT](@ref pm_distPois::LAMBDA_LIMIT) \f$\leq \lambda\f$)
    !>  \param[in]      sqrtLambda      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                      of the same type and kind as `lambda`, representing `sqrt(lambda)`.<br>
    !>                                      This argument is required to ensure fast computation of many random values from the same distribution.<br>
    !>                                      (**optional**. It must be present **if and only if** [LAMBDA_LIMIT](@ref pm_distPois::LAMBDA_LIMIT) \f$\leq \lambda\f$)
    !>  \param[in]      expNegLambda    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                      of the same type and kind as `lambda`, representing `exp(-lambda)`.<br>
    !>                                      (**optional**. It must be present **if and only if** \f$0 < \lambda <\f$ [LAMBDA_LIMIT](@ref pm_distPois::LAMBDA_LIMIT))
    !>
    !>  \interface{setPoisRand}
    !>  \code{.F90}
    !>
    !>      use pm_distPois, only: setPoisRand
    !>
    !>      ! when `LAMBDA_LIMIT <= lambda`
    !>
    !>      call setPoisRand(rand, lambda, logLambda, sqrtLambda)
    !>      call setPoisRand(rand(..), lambda, logLambda, sqrtLambda)
    !>      call setPoisRand(rng, rand, lambda, logLambda, sqrtLambda)
    !>      call setPoisRand(rng, rand(:), lambda, logLambda, sqrtLambda)
    !>
    !>      ! when `0. < lambda < LAMBDA_LIMIT`
    !>
    !>      call setPoisRand(rand, expNegLambda)
    !>      call setPoisRand(rand(..), expNegLambda)
    !>      call setPoisRand(rng, rand, expNegLambda)
    !>      call setPoisRand(rng, rand(:), expNegLambda)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `expNegLambda < 1.` must hold for the corresponding input arguments.<br>
    !>  The condition `exp(-LAMBDA_LIMIT) < expNegLambda` must hold for the corresponding input arguments.<br>
    !>  The condition `logLambda == log(lambda)` must hold for the corresponding input arguments.<br>
    !>  The condition `sqrtLambda == sqrt(lambda)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPoisRand](@ref pm_distPois::getPoisRand)<br>
    !>
    !>  \example{setPoisRand}
    !>  \include{lineno} example/pm_distPois/setPoisRand/main.F90
    !>  \compilef{setPoisRand}
    !>  \output{setPoisRand}
    !>  \include{lineno} example/pm_distPois/setPoisRand/main.out.F90
    !>  \postproc{setPoisRand}
    !>  \include{lineno} example/pm_distPois/setPoisRand/main.py
    !>  \vis{setPoisRand}
    !>  \image html pm_distPois/setPoisRand/setPoisRand.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distPois](@ref test_pm_distPois)
    !>
    !>  \final{setPoisRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setPoisRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setPoisRandExpRNGD_D0_RK5(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setPoisRandExpRNGD_D0_RK4(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setPoisRandExpRNGD_D0_RK3(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setPoisRandExpRNGD_D0_RK2(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setPoisRandExpRNGD_D0_RK1(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setPoisRandExpRNGF_D0_RK5(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setPoisRandExpRNGF_D0_RK4(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setPoisRandExpRNGF_D0_RK3(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setPoisRandExpRNGF_D0_RK2(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setPoisRandExpRNGF_D0_RK1(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPoisRandExpRNGX_D0_RK5(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPoisRandExpRNGX_D0_RK4(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPoisRandExpRNGX_D0_RK3(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPoisRandExpRNGX_D0_RK2(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPoisRandExpRNGX_D0_RK1(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPoisRandExpRNGD_D1_RK5(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPoisRandExpRNGD_D1_RK4(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPoisRandExpRNGD_D1_RK3(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPoisRandExpRNGD_D1_RK2(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPoisRandExpRNGD_D1_RK1(rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPoisRandExpRNGF_D1_RK5(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPoisRandExpRNGF_D1_RK4(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPoisRandExpRNGF_D1_RK3(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPoisRandExpRNGF_D1_RK2(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPoisRandExpRNGF_D1_RK1(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGF_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPoisRandExpRNGX_D1_RK5(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPoisRandExpRNGX_D1_RK4(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPoisRandExpRNGX_D1_RK3(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPoisRandExpRNGX_D1_RK2(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPoisRandExpRNGX_D1_RK1(rng, rand, expNegLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandExpRNGX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: expNegLambda
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
    impure elemental module subroutine setPoisRandRejRNGD_D0_RK5(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setPoisRandRejRNGD_D0_RK4(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setPoisRandRejRNGD_D0_RK3(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setPoisRandRejRNGD_D0_RK2(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setPoisRandRejRNGD_D0_RK1(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setPoisRandRejRNGF_D0_RK5(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setPoisRandRejRNGF_D0_RK4(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setPoisRandRejRNGF_D0_RK3(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setPoisRandRejRNGF_D0_RK2(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setPoisRandRejRNGF_D0_RK1(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setPoisRandRejRNGX_D0_RK5(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setPoisRandRejRNGX_D0_RK4(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setPoisRandRejRNGX_D0_RK3(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setPoisRandRejRNGX_D0_RK2(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setPoisRandRejRNGX_D0_RK1(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPoisRandRejRNGD_D1_RK5(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPoisRandRejRNGD_D1_RK4(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPoisRandRejRNGD_D1_RK3(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPoisRandRejRNGD_D1_RK2(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPoisRandRejRNGD_D1_RK1(rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPoisRandRejRNGF_D1_RK5(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPoisRandRejRNGF_D1_RK4(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPoisRandRejRNGF_D1_RK3(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPoisRandRejRNGF_D1_RK2(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPoisRandRejRNGF_D1_RK1(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGF_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPoisRandRejRNGX_D1_RK5(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPoisRandRejRNGX_D1_RK4(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPoisRandRejRNGX_D1_RK3(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPoisRandRejRNGX_D1_RK2(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPoisRandRejRNGX_D1_RK1(rng, rand, lambda, logLambda, sqrtLambda)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPoisRandRejRNGX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKG)               , intent(in)                    :: lambda, logLambda, sqrtLambda
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distPois
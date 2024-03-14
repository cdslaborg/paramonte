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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Cyclic Geometric distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Cyclic Geometric distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Mass Function (**PMF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the random number generation from the distribution (**RNG**)
!>      <li>    the **Inverse Cumulative Distribution Function (ICDF)** or the **Quantile Function**
!>  </ol>
!>
!>  \details
!>  The Cyclic Geometric distribution is defined similar to the [Geometric distribution](@ref pm_distGeom),
!>  except for the fact that there is an upper limit to the number of Bernoulli trials in the experiment.<br>
!>  Once the upper limit is reached without any success, the experiment recycles to the first Bernoulli trial
!>  and the process repeats until the first success occurs.<br>
!>  The Cyclic Geometric distribution appears in the computation of the workload of the different
!>  processes in a parallel application, for example, the ParaMonte sampler parallel simulations.<br>
!>
!>  **Probability Mass Function (PMF)**<br>
!>
!>  If the probability of success on each trial is \f$\ms{probSuccess}\f$,
!>  then the probability that the \f$\ms{stepSuccess}^{th}\f$ trial is the first success in a cyclic trial
!>  set with a cycle \f$\ms{period}\f$ can be written in terms of the **PMF** of the Geometric distribution as,
!>  \f{eqnarray}{
!>      \large
!>      \pi_{\mathcal{CG}} (X = \ms{stepSuccess} ~|~ \ms{probSuccess}, \ms{period})
!>      &=& \sum_{i = 0}^{+\infty} ~ \pi_{\mathcal{G}} (X = i \times \ms{period} + \ms{stepSuccess} ~|~ \ms{probSuccess}) ~, \nonumber \\
!>      &=& \frac{\ms{probSuccess} (1 - \ms{probSuccess})^{\ms{stepSuccess} - 1}}{1 - (1 - \ms{probSuccess})^{\ms{period}}} ~,
!>  \f}
!>  where,
!>  <ol>
!>      <li> \f$\pi_{\mathcal{CG}} (\cdot)\f$ refers to the [Cyclic Geometric distribution](@ref pm_distGeomCyclic) and,
!>      <li> \f$\pi_{\mathcal{G}} (\cdot)\f$ refers to the [Geometric distribution](@ref pm_distGeom),
!>  </ol>
!>
!>  **Cumulative Distribution Function (CDF)**<br>
!>
!>  The **CDF** of the distribution can be computed as a finite Geometric series,
!>  \f{eqnarray}{
!>      \ms{CDF}(X = \ms{stepSuccess} ~|~ \ms{probSuccess}, \ms{period})
!>      &=& \sum_{i = 0}^{\ms{period}} ~ \pi_{\mathcal{CG}} (X = i ~|~ \ms{probSuccess}, \ms{period}) ~, \nonumber \\
!>      &=& \frac{1 - (1 - \ms{probSuccess})^{\ms{stepSuccess}}}{1 - (1 - \ms{probSuccess})^{\ms{period}}} ~,
!>  \f}
!>
!>  See [Amir Shahmoradi, Fatemeh Bagheri (2020). ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations.](https://arxiv.org/abs/2008.09589)
!>  for details of the derivation of the above PMF.<br>
!>
!>  \see
!>  [pm_distGeom](@ref pm_distGeom)<br>
!>  [Amir Shahmoradi, Fatemeh Bagheri (2020). ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations.](https://arxiv.org/abs/2008.09589)<br>
!>
!>  \test
!>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distGeomCyclic

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distGeomCyclic"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Cyclic Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[0, \ms{period}]\f$.
    !>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  <ol>
    !>                                      <li>    type `real` of the same kind as that of `probSuccess`,<br>
    !>                                      <li>    type `integer` of default kind \IK,<br>
    !>                                  </ol>
    !>                                  containing the values at which the PMF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>                                  As such, it must never exceed the value of the input argument `period`.<br>
    !>  \param[in]  probSuccess     :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the **probability of success** at each Bernoulli step.<br>
    !>  \param[in]  period          :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the period of the Cyclic Geometric distribution.<br>
    !>                                  The `period` represents the maximum number of Bernoulli trials in the experiment.<br>
    !>
    !>  \return
    !>  `logPMF`                    :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `probSuccess`, containing the natural logarithm of the PMF of
    !>                                  the distribution at the specified point `stepSuccess`.<br>
    !>
    !>  \interface{getGeomCyclicLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: getGeomCyclicLogPMF
    !>
    !>      logPMF = getGeomCyclicLogPMF(stepSuccess, probSuccess, period)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < probSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= stepSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `stepSuccess <= period` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomCyclicLogCDF](@ref pm_distGeomCyclic::getGeomCyclicLogCDF)<br>
    !>  [setGeomCyclicLogCDF](@ref pm_distGeomCyclic::setGeomCyclicLogCDF)<br>
    !>  [getGeomCyclicLogPMF](@ref pm_distGeomCyclic::getGeomCyclicLogPMF)<br>
    !>  [setGeomCyclicLogPMF](@ref pm_distGeomCyclic::setGeomCyclicLogPMF)<br>
    !>
    !>  \example{getGeomCyclicLogPMF}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicLogPMF/main.F90
    !>  \compilef{getGeomCyclicLogPMF}
    !>  \output{getGeomCyclicLogPMF}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicLogPMF/main.out.F90
    !>  \postproc{getGeomCyclicLogPMF}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicLogPMF/main.py
    !>  \vis{getGeomCyclicLogPMF}
    !>  \image html pm_distGeomCyclic/getGeomCyclicLogPMF/getGeomCyclicLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \finmain{getGeomCyclicLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getGeomCyclicLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_IK_RK5(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_IK_RK4(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_IK_RK3(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_IK_RK2(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_IK_RK1(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_RK_RK5(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_RK_RK4(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_RK_RK3(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_RK_RK2(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGeomCyclicLogPMF_D0_RK_RK1(stepSuccess, probSuccess, period) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogPMF_D0_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Cyclic Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[0, period]\f$.
    !>
    !>  \param[out] logPMF          :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `logProbSuccess`,
    !>                                  containing the natural logarithm of the PMF of the distribution at the specified `stepSuccess`.<br>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  <ol>
    !>                                      <li>    type `real` of the same kind as that of `logProbSuccess`,<br>
    !>                                      <li>    type `integer` of default kind \IK,<br>
    !>                                  </ol>
    !>                                  containing the values at which the PMF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>  \param[in]  logProbSuccess  :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the natural logarithm of the **probability of success** at each Bernoulli step.<br>
    !>  \param[in]  logProbFailure  :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `logProbSuccess`
    !>                                  containing the natural logarithm of the **probability of failure** at each Bernoulli step.<br>
    !>                                  (**optional**, default = `log(get1mexp(logProbSuccess))`)
    !>
    !>  \interface{setGeomCyclicLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: setGeomCyclicLogPMF
    !>
    !>      call setGeomCyclicLogPMF(logPMF, stepSuccess, logProbSuccess, period) ! elemental
    !>      call setGeomCyclicLogPMF(logPMF, stepSuccess, logProbSuccess, logProbFailure, period) ! elemental
    !>
    !>      call setGeomCyclicLogPMF(logPMF(:), stepSuccess(:), logProbSuccess, period)
    !>      call setGeomCyclicLogPMF(logPMF(:), stepSuccess(:), logProbSuccess, logProbFailure, period)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 < [stepSuccess])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([stepSuccess] <= period)` must hold for the corresponding input arguments.<br>
    !>  The condition `logProbSuccess <= 0.` must hold for the corresponding input arguments.<br>
    !>  The condition `logProbFailure <  0.` must hold for the corresponding input arguments.<br>
    !>  The condition `exp(logProbFailure) + exp(logProbSuccess) == 1.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomCyclicLogCDF](@ref pm_distGeomCyclic::getGeomCyclicLogCDF)<br>
    !>  [setGeomCyclicLogCDF](@ref pm_distGeomCyclic::setGeomCyclicLogCDF)<br>
    !>  [getGeomCyclicLogPMF](@ref pm_distGeomCyclic::getGeomCyclicLogPMF)<br>
    !>  [setGeomCyclicLogPMF](@ref pm_distGeomCyclic::setGeomCyclicLogPMF)<br>
    !>
    !>  \example{setGeomCyclicLogPMF}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicLogPMF/main.F90
    !>  \compilef{setGeomCyclicLogPMF}
    !>  \output{setGeomCyclicLogPMF}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicLogPMF/main.out.F90
    !>  \postproc{setGeomCyclicLogPMF}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicLogPMF/main.py
    !>  \vis{setGeomCyclicLogPMF}
    !>  \image html pm_distGeomCyclic/setGeomCyclicLogPMF/setGeomCyclicLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setGeomCyclicLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! IK

    interface setGeomCyclicLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_IK_RK5(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_IK_RK4(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_IK_RK3(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_IK_RK2(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_IK_RK1(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_IK_RK5(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_IK_RK4(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_IK_RK3(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_IK_RK2(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_IK_RK1(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_IK_RK5(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_IK_RK4(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_IK_RK3(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_IK_RK2(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_IK_RK1(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_IK_RK5(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_IK_RK4(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_IK_RK3(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_IK_RK2(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_IK_RK1(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RK

    interface setGeomCyclicLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_RK_RK5(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_RK_RK4(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_RK_RK3(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_RK_RK2(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFDef_D0_RK_RK1(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D0_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_RK_RK5(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_RK_RK4(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_RK_RK3(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_RK_RK2(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogPMFLog_D0_RK_RK1(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D0_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_RK_RK5(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_RK_RK4(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_RK_RK3(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_RK_RK2(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogPMFDef_D1_RK_RK1(logPMF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFDef_D1_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_RK_RK5(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_RK_RK4(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_RK_RK3(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_RK_RK2(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogPMFLog_D1_RK_RK1(logPMF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogPMFLog_D1_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logPMF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Cumulative Distribution Function (CDF) of the
    !>  <b>Cyclic Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[0, period]\f$.
    !>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  <ol>
    !>                                      <li>    type `real` of the same kind as that of `probSuccess`,<br>
    !>                                      <li>    type `integer` of default kind \IK,<br>
    !>                                  </ol>
    !>                                  containing the values at which the CDF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>                                  As such, it must never exceed the value of the input argument `period`.<br>
    !>  \param[in]  probSuccess     :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the **probability of success** at each Bernoulli step.<br>
    !>  \param[in]  period          :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the period of the Cyclic Geometric distribution.<br>
    !>                                  The `period` represents the maximum number of Bernoulli trials in the experiment.<br>
    !>
    !>  \return
    !>  `logCDF`                    :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `probSuccess`, containing the natural logarithm of the CDF of
    !>                                  the distribution at the specified point `stepSuccess`.<br>
    !>
    !>  \interface{getGeomCyclicLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: getGeomCyclicLogCDF
    !>
    !>      logCDF = getGeomCyclicLogCDF(stepSuccess, probSuccess, period)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < probSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= stepSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `stepSuccess <= period` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomCyclicLogCDF](@ref pm_distGeomCyclic::getGeomCyclicLogCDF)<br>
    !>  [setGeomCyclicLogCDF](@ref pm_distGeomCyclic::setGeomCyclicLogCDF)<br>
    !>  [getGeomCyclicLogPMF](@ref pm_distGeomCyclic::getGeomCyclicLogPMF)<br>
    !>  [setGeomCyclicLogPMF](@ref pm_distGeomCyclic::setGeomCyclicLogPMF)<br>
    !>
    !>  \example{getGeomCyclicLogCDF}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicLogCDF/main.F90
    !>  \compilef{getGeomCyclicLogCDF}
    !>  \output{getGeomCyclicLogCDF}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicLogCDF/main.out.F90
    !>  \postproc{getGeomCyclicLogCDF}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicLogCDF/main.py
    !>  \vis{getGeomCyclicLogCDF}
    !>  \image html pm_distGeomCyclic/getGeomCyclicLogCDF/getGeomCyclicLogCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \finmain{getGeomCyclicLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getGeomCyclicLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_IK_RK5(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_IK_RK4(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_IK_RK3(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_IK_RK2(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_IK_RK1(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_RK_RK5(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_RK_RK4(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_RK_RK3(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_RK_RK2(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGeomCyclicLogCDF_D0_RK_RK1(stepSuccess, probSuccess, period) result(logCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicLogCDF_D0_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Cumulative Distribution Function (CDF) of the
    !>  <b>Cyclic Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[0, period]\f$.
    !>
    !>  \param[out] logCDF          :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `logProbSuccess`,
    !>                                  containing the natural logarithm of the CDF of the distribution at the specified `stepSuccess`.<br>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  <ol>
    !>                                      <li>    type `real` of the same kind as that of `logProbSuccess`,<br>
    !>                                      <li>    type `integer` of default kind \IK,<br>
    !>                                  </ol>
    !>                                  containing the values at which the CDF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>  \param[in]  logProbSuccess  :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the natural logarithm of the **probability of success** at each Bernoulli step.<br>
    !>  \param[in]  logProbFailure  :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `logProbSuccess`
    !>                                  containing the natural logarithm of the **probability of failure** at each Bernoulli step.<br>
    !>                                  (**optional**, default = `log(get1mexp(logProbSuccess))`)
    !>  \param[in]  period          :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the period of the Cyclic Geometric distribution.<br>
    !>                                  The `period` represents the maximum number of Bernoulli trials in the experiment.<br>
    !>
    !>  \interface{setGeomCyclicLogCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: setGeomCyclicLogCDF
    !>
    !>      call setGeomCyclicLogCDF(logCDF, stepSuccess, logProbSuccess, period) ! elemental
    !>      call setGeomCyclicLogCDF(logCDF, stepSuccess, logProbSuccess, logProbFailure, period) ! elemental
    !>
    !>      call setGeomCyclicLogCDF(logCDF(:), stepSuccess(:), logProbSuccess, period)
    !>      call setGeomCyclicLogCDF(logCDF(:), stepSuccess(:), logProbSuccess, logProbFailure, period)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 < [stepSuccess])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([stepSuccess] <= period)` must hold for the corresponding input arguments.<br>
    !>  The condition `logProbSuccess <= 0.` must hold for the corresponding input arguments.<br>
    !>  The condition `logProbFailure <  0.` must hold for the corresponding input arguments.<br>
    !>  The condition `exp(logProbFailure) + exp(logProbSuccess) == 1.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomCyclicLogCDF](@ref pm_distGeomCyclic::getGeomCyclicLogCDF)<br>
    !>  [setGeomCyclicLogCDF](@ref pm_distGeomCyclic::setGeomCyclicLogCDF)<br>
    !>  [getGeomCyclicLogPMF](@ref pm_distGeomCyclic::getGeomCyclicLogPMF)<br>
    !>  [setGeomCyclicLogPMF](@ref pm_distGeomCyclic::setGeomCyclicLogPMF)<br>
    !>
    !>  \example{setGeomCyclicLogCDF}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicLogCDF/main.F90
    !>  \compilef{setGeomCyclicLogCDF}
    !>  \output{setGeomCyclicLogCDF}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicLogCDF/main.out.F90
    !>  \postproc{setGeomCyclicLogCDF}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicLogCDF/main.py
    !>  \vis{setGeomCyclicLogCDF}
    !>  \image html pm_distGeomCyclic/setGeomCyclicLogCDF/setGeomCyclicLogCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setGeomCyclicLogCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! IK

    interface setGeomCyclicLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_IK_RK5(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_IK_RK4(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_IK_RK3(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_IK_RK2(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_IK_RK1(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_IK_RK5(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_IK_RK4(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_IK_RK3(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_IK_RK2(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_IK_RK1(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: stepSuccess, period
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_IK_RK5(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_IK_RK4(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_IK_RK3(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_IK_RK2(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_IK_RK1(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_IK_RK5(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_IK_RK4(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_IK_RK3(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_IK_RK2(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_IK_RK1(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RK

    interface setGeomCyclicLogCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_RK_RK5(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_RK_RK4(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_RK_RK3(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_RK_RK2(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFDef_D0_RK_RK1(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D0_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_RK_RK5(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_RK_RK4(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_RK_RK3(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_RK_RK2(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicLogCDFLog_D0_RK_RK1(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D0_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_RK_RK5(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_RK_RK4(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_RK_RK3(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_RK_RK2(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogCDFDef_D1_RK_RK1(logCDF, stepSuccess, logProbSuccess, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFDef_D1_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_RK_RK5(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_RK_RK4(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_RK_RK3(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_RK_RK2(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicLogCDFLog_D1_RK_RK1(logCDF, stepSuccess, logProbSuccess, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicLogCDFLog_D1_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)   , contiguous    :: logCDF(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank of) random value(s) from the Cyclic Geometric distribution.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distGeomCyclic](@ref pm_distGeomCyclic) for more details.<br>
    !>
    !>  \param[in]  probSuccess :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              representing the probability of success parameter of the distribution.<br>
    !>  \param[in]  period      :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the period of the Cyclic Geometric distribution.<br>
    !>                              The `period` represents the maximum number of Bernoulli trials in the experiment.<br>
    !>
    !>  \return
    !>  `rand`                  :   The output positive scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `probSuccess`, containing random value(s) from the distribution.<br>
    !>
    !>  \interface{getGeomCyclicRand}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: getGeomCyclicRand
    !>
    !>      rand = getGeomCyclicRand(probSuccess, period)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < period` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < probSuccess .and. probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGeomCyclicRand](@ref pm_distGeomCyclic::setGeomCyclicRand)<br>
    !>
    !>  \example{getGeomCyclicRand}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicRand/main.F90
    !>  \compilef{getGeomCyclicRand}
    !>  \output{getGeomCyclicRand}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicRand/main.out.F90
    !>  \postproc{getGeomCyclicRand}
    !>  \include{lineno} example/pm_distGeomCyclic/getGeomCyclicRand/main.py
    !>  \vis{getGeomCyclicRand}
    !>  \image html pm_distGeomCyclic/getGeomCyclicRand/getGeomCyclicRand.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \finmain{getGeomCyclicRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGeomCyclicRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGeomCyclicRand_RK5(probSuccess, period) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicRand_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                                 :: rand
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGeomCyclicRand_RK4(probSuccess, period) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicRand_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                                 :: rand
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGeomCyclicRand_RK3(probSuccess, period) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicRand_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                                 :: rand
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGeomCyclicRand_RK2(probSuccess, period) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicRand_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                                 :: rand
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGeomCyclicRand_RK1(probSuccess, period) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCyclicRand_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                                 :: rand
        integer(IK) , intent(in)                    :: period
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Cyclic Geometric distribution.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distGeomCyclic](@ref pm_distGeomCyclic) for more details.<br>
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
    !>  \param[in]      logProbFailure  :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                      <ol>
    !>                                          <li>    type `real` of kind \RKALL,
    !>                                      </ol>
    !>                                      representing the natural logarithm of the probability of failure(s) `log(1 - probSuccess, period)` of the distribution(s).<br>
    !>  \param[in]      period          :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                      type `integer` of default kind \IK, containing the period of the Cyclic Geometric distribution.<br>
    !>                                      The `period` represents the maximum number of Bernoulli trials in the experiment.<br>
    !>
    !>  \interface{setGeomCyclicRand}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: setGeomCyclicRand
    !>      use pm_kind, only: IK
    !>      integer(IK) :: rand
    !>
    !>      call setGeomCyclicRand(rand, logProbFailure, period) ! elemental
    !>      call setGeomCyclicRand(rng, rand, logProbFailure, period)
    !>      call setGeomCyclicRand(rng, rand(:), logProbFailure(:), period)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < period` must hold for the corresponding input arguments.<br>
    !>  The condition `logProbFailure < 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomCyclicRand](@ref pm_distGeomCyclic::getGeomCyclicRand)<br>
    !>
    !>  \example{setGeomCyclicRand}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicRand/main.F90
    !>  \compilef{setGeomCyclicRand}
    !>  \output{setGeomCyclicRand}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicRand/main.out.F90
    !>  \postproc{setGeomCyclicRand}
    !>  \include{lineno} example/pm_distGeomCyclic/setGeomCyclicRand/main.py
    !>  \vis{setGeomCyclicRand}
    !>  \image html pm_distGeomCyclic/setGeomCyclicRand/setGeomCyclicRand.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \finmain{setGeomCyclicRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGeomCyclicRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGD_D0_RK5(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGD_D0_RK4(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGD_D0_RK3(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGD_D0_RK2(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGD_D0_RK1(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGF_D0_RK5(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGF_D0_RK4(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGF_D0_RK3(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGF_D0_RK2(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGeomCyclicRandRNGF_D0_RK1(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCyclicRandRNGX_D0_RK5(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCyclicRandRNGX_D0_RK4(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCyclicRandRNGX_D0_RK3(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCyclicRandRNGX_D0_RK2(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCyclicRandRNGX_D0_RK1(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setGeomCyclicRandRNGD_D1_RK5(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setGeomCyclicRandRNGD_D1_RK4(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setGeomCyclicRandRNGD_D1_RK3(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setGeomCyclicRandRNGD_D1_RK2(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setGeomCyclicRandRNGD_D1_RK1(rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGD_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setGeomCyclicRandRNGF_D1_RK5(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setGeomCyclicRandRNGF_D1_RK4(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setGeomCyclicRandRNGF_D1_RK3(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setGeomCyclicRandRNGF_D1_RK2(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setGeomCyclicRandRNGF_D1_RK1(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGF_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomCyclicRandRNGX_D1_RK5(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomCyclicRandRNGX_D1_RK4(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomCyclicRandRNGX_D1_RK3(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomCyclicRandRNGX_D1_RK2(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomCyclicRandRNGX_D1_RK1(rng, rand, logProbFailure, period)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCyclicRandRNGX_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        integer(IK)             , intent(in)                    :: period
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the parameters of a least-squares **fit to the histogram**
    !>  representing a Cyclic-Geometric-distributed sample can be successfully inferred, otherwise, return `.false.`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distGeomCyclic](@ref pm_distGeomCyclic) for more details on the Cyclic Geometric distribution.<br>
    !>  Given a set of Cyclic Bernoulli trial success steps \f$\ms{stepSuccess}\f$ and the corresponding number of successes \f$\ms{freqSuccess}\f$
    !>  at these steps, with the supplied Cyclic Bernoulli trial \f$\ms{period}\f$, the procedures of this generic interface return
    !>  the best-fit parameters \f$(\ms{probSuccess}, \ms{normFac})\f$ of the following curve,
    !>  \f{eqnarray}{
    !>      \large
    !>      \ms{freqSuccess}
    !>      &=& \ms{normFac} \times \pi_{\mathcal{CG}} (X = \ms{stepSuccess} ~|~ \ms{probSuccess}, \ms{period}) ~, \nonumber \\
    !>      &=& \ms{normFac} \times \frac{\ms{probSuccess} (1 - \ms{probSuccess})^{\ms{stepSuccess} - 1}}{1 - (1 - \ms{probSuccess})^{\ms{period}}} ~,
    !>  \f}
    !>  where the probability of success on each trial is \f$\ms{probSuccess}\f$,
    !>  and \f$\pi_{\mathcal{CG}}\f$ is the [probability mass function](@ref pm_distGeomCyclic) of the \f$\ms{stepSuccess}\f$th
    !>  trial being the first success in a cyclic trial set with the cycle \f$\ms{period}\f$.<br>
    !>
    !>  \note
    !>  This fitting appears frequently in Fork-Join parallel MCMC simulations of the ParaMonte library.<br>
    !>  See [Amir Shahmoradi, Fatemeh Bagheri (2020). ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations.](https://arxiv.org/abs/2008.09589)
    !>  for details and discussion.<br>
    !>
    !>  \param[in]  stepSuccess     :   The input `contiguous` positive vector of<br>
    !>                                  <ol>
    !>                                      <li>    type `real` of the same kind as that of `probSuccess`,<br>
    !>                                      <li>    type `integer` of default kind \IK,<br>
    !>                                  </ol>
    !>                                  containing the steps (\f$x\f$) at which contribution frequencies represented by input vector `freqSuccess` have occurred.<br>
    !>  \param[in]  freqSuccess     :   The input `contiguous` positive vector of the same type, kind, and size as the input `stepSuccess`,
    !>                                  containing the contribution frequencies at the corresponding `stepSuccess` in the histogram.<br>
    !>  \param[in]  period          :   The input positive scalar of type `integer` of default kind \IK, containing the period of the Cyclic Geometric distribution.<br>
    !>                                  The `period` represents the maximum number of Bernoulli trials in the experiment.<br>
    !>                                  Note that the condition `maxval(stepSuccess) <= period` must hold and the sizes
    !>                                  of `stepSuccess` and `freqSuccess` vectors must be smaller than the specified `period`.<br>
    !>  \param[out] probSuccess     :   The output positive scalar of<br>
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  representing the best-fit probability-of-success parameter of the Cyclic Geometric distribution inferred from the fitting.<br>
    !>  \param[in]  normFac         :   The output positive scalar of the same type and kind as the output `probSuccess`, containing the
    !>                                  best-fit normalization constant of the histogram represented by the input `stepSuccess` and freqSuccess`.<br>
    !>
    !>  \return
    !>  `failed`                    :   The output scalar of type `logical` of default kind \LK that is `.true.`
    !>                                  if and only if the best-fit output parameters are successfully inferred.<br>
    !>                                  The algorithm can fail only if the optimizer fails to converge to set of best-fit parameters.<br>
    !>
    !>  \interface{isFailedGeomCyclicFit}
    !>  \code{.F90}
    !>
    !>      use pm_distGeomCyclic, only: isFailedGeomCyclicFit
    !>
    !>      failed = isFailedGeomCyclicFit(stepSuccess(1:period), freqSuccess(1:period), period, probSuccess, normFac)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 < stepSuccess)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 < size(freqSuccess))` must hold for the corresponding input arguments.<br>
    !>  The condition `size(stepSuccess) == size(freqSuccess)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(stepSuccess <= size(stepSuccess))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(stepSuccess <= period)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setGeomCyclicRand](@ref pm_distGeomCyclic::setGeomCyclicRand)<br>
    !>  [setGeomCyclicLogPMF](@ref pm_distGeomCyclic::setGeomCyclicLogPMF)<br>
    !>  [setGeomCyclicLogCDF](@ref pm_distGeomCyclic::setGeomCyclicLogCDF)<br>
    !>
    !>  \example{isFailedGeomCyclicFit}
    !>  \include{lineno} example/pm_distGeomCyclic/isFailedGeomCyclicFit/main.F90
    !>  \compilef{isFailedGeomCyclicFit}
    !>  \output{isFailedGeomCyclicFit}
    !>  \include{lineno} example/pm_distGeomCyclic/isFailedGeomCyclicFit/main.out.F90
    !>  \postproc{isFailedGeomCyclicFit}
    !>  \include{lineno} example/pm_distGeomCyclic/isFailedGeomCyclicFit/main.py
    !>  \vis{isFailedGeomCyclicFit}
    !>  \image html pm_distGeomCyclic/isFailedGeomCyclicFit/isFailedGeomCyclicFit.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeomCyclic](@ref test_pm_distGeomCyclic)
    !>
    !>  \finmain{isFailedGeomCyclicFit}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    !  \param[in]  freqSuccess :   The input `contiguous` non-negative vector of the same type, kind, and size as the input `stepSuccess`,
    !                              containing the contribution frequencies at the corresponding `stepSuccess` in the histogram.<br>
    interface isFailedGeomCyclicFit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function isFailedGeomCyclicFit_IK_RK5(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_IK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK4_ENABLED
    impure module function isFailedGeomCyclicFit_IK_RK4(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_IK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK3_ENABLED
    impure module function isFailedGeomCyclicFit_IK_RK3(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_IK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK2_ENABLED
    impure module function isFailedGeomCyclicFit_IK_RK2(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_IK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK1_ENABLED
    impure module function isFailedGeomCyclicFit_IK_RK1(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_IK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: stepSuccess(:)
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function isFailedGeomCyclicFit_RK_RK5(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_RK_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK4_ENABLED
    impure module function isFailedGeomCyclicFit_RK_RK4(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_RK_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK3_ENABLED
    impure module function isFailedGeomCyclicFit_RK_RK3(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_RK_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK2_ENABLED
    impure module function isFailedGeomCyclicFit_RK_RK2(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_RK_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

#if RK1_ENABLED
    impure module function isFailedGeomCyclicFit_RK_RK1(stepSuccess, freqSuccess, period, probSuccess, normFac) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGeomCyclicFit_RK_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: period
        integer(IK) , intent(in)    , contiguous    :: freqSuccess(:)
        real(RKC)   , intent(in)    , contiguous    :: stepSuccess(:)
        real(RKC)   , intent(out)                   :: probSuccess, normFac
        logical(LK)                                 :: failed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distGeomCyclic
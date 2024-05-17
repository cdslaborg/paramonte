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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Geometric distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Geometric distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Mass Function (**PMF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **Geometric distribution** is either one of two discrete probability distributions:
!>  <ol>
!>      <li>    The probability distribution of the number \f$X\f$ of [Bernoulli trials](@ref pm_distBern) needed to get one success,
!>              supported on the set \f$\{1,2,3,\ldots\}\f$.
!>      <li>    The probability distribution of the number \f$Y = X − 1\f$ of failures before the first success,
!>              supported on the set \f$\{0,1,2,\ldots \}\f$.
!>  </ol>
!>  Which of these is called the Geometric distribution is a matter of convention and convenience.<br>
!>  These two different Geometric distributions should not be confused with each other.<br>
!>  Frequently, the name **shifted Geometric distribution** is adopted for the former (distribution of the number \f$X\f$).<br>
!>  However, to avoid ambiguity, it is considered wise to indicate which is intended, by mentioning the support explicitly.<br>
!>  <b>The generic interfaces of this module return the probability mass function of the number \f$X\f$, i.e.,
!>  the probability of the \f$X\f$ number of Bernoulli trials needed to get one success</b>.<br>
!>  The Geometric distribution gives the probability that the first occurrence of success requires \f$k\f$ independent trials,
!>  each with success probability \f$p\f$.<br>
!>  If the probability of success on each trial is \f$p\f$, then the probability that the \f$k\f$th trial is the first success is
!>  \f{equation}{
!>      \pi(X = k | p) = (1 - p)^{k - 1} p ~,~ k = 1, 2, 3, 4, \ldots ~.
!>  \f}
!>  The above form of the Geometric distribution is used for modeling the number of trials up to and including the first success.<br>
!>  <b>This is the form implemented in this module</b>.<br>
!>  By contrast, the following form of the geometric distribution is used for modeling the number of failures until the first success:
!>  \f{equation}{
!>      \pi(Y = k) = \pi(X = k + 1) = (1 - p)^k p ~,~ k = 0, 1, 2, 3, \ldots ~.
!>  \f}
!>  In either case, the sequence of probabilities is a [geometric sequence](https://en.wikipedia.org/wiki/Geometric_progression).<br>
!>  **Cumulative Distribution Function (CDF)**<br>
!>
!>  The **CDF** of the Geometric distribution with parameter \f$\ms{probSuccess}\f$ (representing the probability of success) is defined as,
!>  \f{equation}{
!>      \ms{CDF}(\ms{stepSuccess} ~|~ \ms{probSuccess}) = 1 - (1 - \ms{probSuccess})^{\ms{stepSuccess}} ~,~ 0 < \ms{stepSuccess} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$\ms{stepSuccess}\f$ is the number of Bernoulli trial steps up to and including the first success.<br>
!>  </ol>
!>
!>  **Random Number Generation (RNG)**<br>
!>
!>  The [exponential distribution](@ref pm_distExp) is the continuous analogue of the [geometric distribution](@ref pm_distGeom).<br>
!>  If \f$X\f$ is an exponentially distributed random variable with parameter \f$\lambda\f$, then
!>  \f{equation}{
!>      Y = \lfloor X \rfloor ~,
!>  \f}
!>  where \f$\lfloor\quad\rfloor\f$ is the `floor()` (or greatest integer) function, is a geometrically distributed random variable
!>  with parameter \f$\ms{probSuccess} = 1 − e^{−\lambda}\f$ (thus \f$\lambda = −ln(1 − \ms{probSuccess})\f$) and taking values in the set \f$\{0, 1, 2, \ldots\}\f$.<br>
!>  This can be used to generate geometrically distributed pseudorandom numbers.<br>
!>  If \f$U\f$ is uniformly distributed in \f$(0,1]\f$, then \f$1 + \lfloor\ln(U)/\ln(1 - \ms{probSuccess})\rfloor\f$ is geometrically distributed with parameter \f$\ms{probSuccess}\f$.<br>
!>
!>  \see
!>  [pm_distGeomCyclic](@ref pm_distGeomCyclic)<br>
!>  [pm_distLogUnif](@ref pm_distLogUnif)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distBern](@ref pm_distBern)<br>
!>  [pm_distExp](@ref pm_distExp)<br>
!>  [Geometric distribution](https://en.wikipedia.org/wiki/Geometric_distribution)<br>
!>
!>  \test
!>  [test_pm_distGeom](@ref test_pm_distGeom)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distGeom

    use pm_kind, only: SK, IK
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distGeom"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Geometric
    !>  as defined in the description of [pm_distGeom](@ref pm_distGeom).
    !>
    !>  \details
    !>  See the documentation of [pm_distGeom](@ref pm_distGeom) for the definition of the Geometric distribution.
    !>
    !>  \interface{distGeom_type}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: distGeom_type
    !>      type(distGeom_type) :: distGeom
    !>
    !>      distGeom = distGeom_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \finmain{distGeom_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distGeom_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[0, +\infty)\f$.
    !>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the values at which the PMF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>  \param[in]  probSuccess     :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the **probability of success** at each Bernoulli step.<br>
    !>
    !>  \return
    !>  `logPMF`                    :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `probSuccess`, containing the natural logarithm of the PMF of
    !>                                  the distribution at the specified point `stepSuccess`.<br>
    !>
    !>  \interface{getGeomLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: getGeomLogPMF
    !>
    !>      logPMF = getGeomLogPMF(stepSuccess, probSuccess)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= stepSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `0. < probSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGeomLogPMF](@ref pm_distGeom::setGeomLogPMF)<br>
    !>
    !>  \example{getGeomLogPMF}
    !>  \include{lineno} example/pm_distGeom/getGeomLogPMF/main.F90
    !>  \compilef{getGeomLogPMF}
    !>  \output{getGeomLogPMF}
    !>  \include{lineno} example/pm_distGeom/getGeomLogPMF/main.out.F90
    !>  \postproc{getGeomLogPMF}
    !>  \include{lineno} example/pm_distGeom/getGeomLogPMF/main.py
    !>  \vis{getGeomLogPMF}
    !>  \image html pm_distGeom/getGeomLogPMF/getGeomLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \finmain{getGeomLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGeomLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGeomLogPMF_RK5(stepSuccess, probSuccess) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomLogPMF_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGeomLogPMF_RK4(stepSuccess, probSuccess) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomLogPMF_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGeomLogPMF_RK3(stepSuccess, probSuccess) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomLogPMF_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGeomLogPMF_RK2(stepSuccess, probSuccess) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomLogPMF_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGeomLogPMF_RK1(stepSuccess, probSuccess) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomLogPMF_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: logPMF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[0, +\infty)\f$.
    !>
    !>  \param[out] logPMF          :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `logProbSuccess`,
    !>                                  containing the natural logarithm of the PMF of the distribution at the specified `stepSuccess`.<br>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the values at which the PMF must be computed.<br>
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
    !>  \interface{setGeomLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: setGeomLogPMF
    !>
    !>      call setGeomLogPMF(logPMF, stepSuccess, logProbSuccess)
    !>      call setGeomLogPMF(logPMF, stepSuccess, logProbSuccess, logProbFailure)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < stepSuccess` must hold for the corresponding input arguments.<br>
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
    !>  [getGeomLogPMF](@ref pm_distGeom::getGeomLogPMF)<br>
    !>
    !>  \example{setGeomLogPMF}
    !>  \include{lineno} example/pm_distGeom/setGeomLogPMF/main.F90
    !>  \compilef{setGeomLogPMF}
    !>  \output{setGeomLogPMF}
    !>  \include{lineno} example/pm_distGeom/setGeomLogPMF/main.out.F90
    !>  \postproc{setGeomLogPMF}
    !>  \include{lineno} example/pm_distGeom/setGeomLogPMF/main.py
    !>  \vis{setGeomLogPMF}
    !>  \image html pm_distGeom/setGeomLogPMF/setGeomLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain{setGeomLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGeomLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomLogPMFDef_RK5(logPMF, stepSuccess, logProbSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFDef_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomLogPMFDef_RK4(logPMF, stepSuccess, logProbSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFDef_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomLogPMFDef_RK3(logPMF, stepSuccess, logProbSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFDef_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomLogPMFDef_RK2(logPMF, stepSuccess, logProbSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFDef_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomLogPMFDef_RK1(logPMF, stepSuccess, logProbSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFDef_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: logProbSuccess
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomLogPMFLog_RK5(logPMF, stepSuccess, logProbSuccess, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFLog_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomLogPMFLog_RK4(logPMF, stepSuccess, logProbSuccess, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFLog_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomLogPMFLog_RK3(logPMF, stepSuccess, logProbSuccess, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFLog_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomLogPMFLog_RK2(logPMF, stepSuccess, logProbSuccess, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFLog_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomLogPMFLog_RK1(logPMF, stepSuccess, logProbSuccess, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomLogPMFLog_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: logProbSuccess, logProbFailure
        real(RKC)   , intent(out)                   :: logPMF
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  <b>Geometric distribution</b> for an input `stepSuccess` within the discrete integer support of the distribution \f$[1, +\infty)\f$.
    !>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the values at which the CDF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>  \param[in]  probSuccess     :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the **probability of success** at each Bernoulli step.<br>
    !>
    !>  \return
    !>  `cdf`                       :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `probSuccess`, containing the natural logarithm of the CDF of
    !>                                  the distribution at the specified point.<br>
    !>
    !>  \interface{getGeomCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: getGeomCDF
    !>
    !>      cdf = getGeomCDF(stepSuccess, probSuccess)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < stepSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < probSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGeomCDF](@ref pm_distGeom::setGeomCDF)<br>
    !>
    !>  \example{getGeomCDF}
    !>  \include{lineno} example/pm_distGeom/getGeomCDF/main.F90
    !>  \compilef{getGeomCDF}
    !>  \output{getGeomCDF}
    !>  \include{lineno} example/pm_distGeom/getGeomCDF/main.out.F90
    !>  \postproc{getGeomCDF}
    !>  \include{lineno} example/pm_distGeom/getGeomCDF/main.py
    !>  \vis{getGeomCDF}
    !>  \image html pm_distGeom/getGeomCDF/getGeomCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \finmain{getGeomCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGeomCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getGeomCDF_RK5(stepSuccess, probSuccess) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getGeomCDF_RK4(stepSuccess, probSuccess) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getGeomCDF_RK3(stepSuccess, probSuccess) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getGeomCDF_RK2(stepSuccess, probSuccess) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getGeomCDF_RK1(stepSuccess, probSuccess) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK) , intent(in)                    :: stepSuccess
        real(RKC)   , intent(in)                    :: probSuccess
        real(RKC)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the <b>Geometric distribution</b>.
    !>
    !>  \param[out] cdf             :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                                  of the same type and kind as `probSuccess`,
    !>                                  containing the natural logarithm of the CDF of the distribution at the specified input value(s).<br>
    !>  \param[in]  stepSuccess     :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                                  type `integer` of default kind \IK, containing the values at which the CDF must be computed.<br>
    !>                                  Note that `stepSuccess` represents the Bernoulli trial step at which the first success occurs.<br>
    !>  \param[in]  probSuccess     :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the **probability of success** at each Bernoulli step.<br>
    !>
    !>  \interface{setGeomCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: setGeomCDF
    !>
    !>      call setGeomCDF(cdf, stepSuccess, probSuccess)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < stepSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < probSuccess` must hold for the corresponding input arguments.<br>
    !>  The condition `probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomCDF](@ref pm_distGeom::getGeomCDF)<br>
    !>
    !>  \example{setGeomCDF}
    !>  \include{lineno} example/pm_distGeom/setGeomCDF/main.F90
    !>  \compilef{setGeomCDF}
    !>  \output{setGeomCDF}
    !>  \include{lineno} example/pm_distGeom/setGeomCDF/main.out.F90
    !>  \postproc{setGeomCDF}
    !>  \include{lineno} example/pm_distGeom/setGeomCDF/main.py
    !>  \vis{setGeomCDF}
    !>  \image html pm_distGeom/setGeomCDF/setGeomCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \finmain{setGeomCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGeomCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomCDF_RK5(cdf, stepSuccess, probSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: probSuccess
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomCDF_RK4(cdf, stepSuccess, probSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: probSuccess
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomCDF_RK3(cdf, stepSuccess, probSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: probSuccess
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomCDF_RK2(cdf, stepSuccess, probSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: probSuccess
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomCDF_RK1(cdf, stepSuccess, probSuccess)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomCDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: probSuccess
        integer(IK) , intent(in)                    :: stepSuccess
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank of) random value(s) from the Geometric distribution.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distGeom](@ref pm_distGeom) for more details.<br>
    !>
    !>  \param[in]  probSuccess :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              representing the probability of success parameter of the distribution.<br>
    !>
    !>  \return
    !>  `rand`                  :   The output positive scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `probSuccess`, containing random value(s) from the distribution.<br>
    !>
    !>  \interface{getGeomRand}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: getGeomRand
    !>
    !>      rand = getGeomRand(probSuccess)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < probSuccess .and. probSuccess <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setGeomRand](@ref pm_distGeom::setGeomRand)<br>
    !>
    !>  \example{getGeomRand}
    !>  \include{lineno} example/pm_distGeom/getGeomRand/main.F90
    !>  \compilef{getGeomRand}
    !>  \output{getGeomRand}
    !>  \include{lineno} example/pm_distGeom/getGeomRand/main.out.F90
    !>  \postproc{getGeomRand}
    !>  \include{lineno} example/pm_distGeom/getGeomRand/main.py
    !>  \vis{getGeomRand}
    !>  \image html pm_distGeom/getGeomRand/getGeomRand.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \finmain{getGeomRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getGeomRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getGeomRand_RK5(probSuccess) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomRand_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                                 :: rand
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getGeomRand_RK4(probSuccess) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomRand_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                                 :: rand
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getGeomRand_RK3(probSuccess) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomRand_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                                 :: rand
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getGeomRand_RK2(probSuccess) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomRand_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                                 :: rand
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getGeomRand_RK1(probSuccess) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeomRand_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                                 :: rand
        real(RKC)   , intent(in)                    :: probSuccess
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Geometric distribution.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distGeom](@ref pm_distGeom) for more details.<br>
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
    !>                                      representing the natural logarithm of the probability of failure(s) `log(1 - probSuccess)` of the distribution(s).<br>
    !>
    !>  \interface{setGeomRand}
    !>  \code{.F90}
    !>
    !>      use pm_distGeom, only: setGeomRand
    !>      use pm_kind, only: IK
    !>      integer(IK) :: rand
    !>
    !>      call setGeomRand(rand, logProbFailure) ! elemental
    !>      call setGeomRand(rng, rand, logProbFailure)
    !>      call setGeomRand(rng, rand(:), logProbFailure(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `logProbFailure < 0.` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getGeomRand](@ref pm_distGeom::getGeomRand)<br>
    !>
    !>  \example{setGeomRand}
    !>  \include{lineno} example/pm_distGeom/setGeomRand/main.F90
    !>  \compilef{setGeomRand}
    !>  \output{setGeomRand}
    !>  \include{lineno} example/pm_distGeom/setGeomRand/main.out.F90
    !>  \postproc{setGeomRand}
    !>  \include{lineno} example/pm_distGeom/setGeomRand/main.py
    !>  \vis{setGeomRand}
    !>  \image html pm_distGeom/setGeomRand/setGeomRand.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distGeom](@ref test_pm_distGeom)
    !>
    !>  \finmain{setGeomRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setGeomRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGeomRandRNGD_D0_RK5(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGeomRandRNGD_D0_RK4(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGeomRandRNGD_D0_RK3(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGeomRandRNGD_D0_RK2(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGeomRandRNGD_D0_RK1(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setGeomRandRNGF_D0_RK5(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setGeomRandRNGF_D0_RK4(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setGeomRandRNGF_D0_RK3(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setGeomRandRNGF_D0_RK2(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setGeomRandRNGF_D0_RK1(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setGeomRandRNGX_D0_RK5(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setGeomRandRNGX_D0_RK4(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setGeomRandRNGX_D0_RK3(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setGeomRandRNGX_D0_RK2(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setGeomRandRNGX_D0_RK1(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setGeomRandRNGD_D1_RK5(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setGeomRandRNGD_D1_RK4(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setGeomRandRNGD_D1_RK3(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setGeomRandRNGD_D1_RK2(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setGeomRandRNGD_D1_RK1(rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGD_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setGeomRandRNGF_D1_RK5(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setGeomRandRNGF_D1_RK4(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setGeomRandRNGF_D1_RK3(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setGeomRandRNGF_D1_RK2(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setGeomRandRNGF_D1_RK1(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGF_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setGeomRandRNGX_D1_RK5(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setGeomRandRNGX_D1_RK4(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setGeomRandRNGX_D1_RK3(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setGeomRandRNGX_D1_RK2(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setGeomRandRNGX_D1_RK1(rng, rand, logProbFailure)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setGeomRandRNGX_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(out)                   :: rand(:)
        real(RKC)               , intent(in)                    :: logProbFailure
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distGeom
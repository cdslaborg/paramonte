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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Binomial distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Binomial distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Mass Function (**PMF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The binomial distribution with parameters \f$n\f$ and \f$p\f$ is the discrete probability distribution of the number of successes
!>  in a sequence of \f$n\f$ independent experiments, each asking a yes–no question, and each with its own Boolean-valued outcome:<br>
!>  success (with probability \f$p\f$) or failure (with probability \f$q = 1 - p\f$).<br>
!>  A single success/failure experiment is also called a Bernoulli trial or Bernoulli experiment,
!>  and a sequence of outcomes is called a Bernoulli process.<br>
!>  For a single trial, i.e., \f$n = 1\f$, the binomial distribution is a [Bernoulli distribution](@ref pm_distBern).<br>
!>  The binomial distribution is the basis for the popular binomial test of statistical significance.<br>
!>  The binomial distribution is frequently used to model the number of successes in a
!>  sample of size \f$n\f$ drawn with replacement from a population of size \f$N\f$.<br>
!>  If the sampling is carried out without replacement, the draws are not independent
!>  and so the resulting distribution is a **hypergeometric distribution**.<br>
!>  However, for \f$N\f$ much larger than \f$n\f$, the binomial distribution
!>  remains a good approximation, and is widely used.<br>
!>
!>  Probability Mass Function
!>  -------------------------
!>
!>  In general, if the random variable \f$X\f$ follows the binomial distribution with parameters \f$\mathbb{N}\f$
!>  and \f$p \in [0, 1]\f$, we write \f$X \sim B(n, p)\f$.<br>
!>  The probability of getting exactly \f$k\f$ successes in \f$n\f$ independent Bernoulli
!>  trials (with the same rate \f$p\f$) is given by the probability mass function:<br>
!>  \f{equation}{
!>      \large
!>      f(k, n, p) = \Pr(X = k) = \binom{n}{k} p^{k}(1-p)^{n-k} ~,
!>  \f}
!>  for \f$k = 0, 1, 2, ..., n\f$, where
!>  \f{equation}{
!>      \large
!>      \binom{n}{k} = \frac{n!}{k!(n-k)!} ~,
!>  \f}
!>  is the **binomial coefficient**, hence the distribution name.<br>
!>  The formula can be understood as follows:<br>
!>  \f$p^k q^{n-k}\f$ is the probability of obtaining the sequence of n Bernoulli trials in which
!>  the first \f$k\f$ trials are *successes* and the remaining last \f$n - k\f$ trials result in *failure*.<br>
!>  Since the trials are independent with probabilities remaining constant between them,
!>  any sequence of \f$n\f$ trials with \f$k\f$ successes and \f$n - k\f$ failures)
!>  has the same probability of being achieved regardless of positions of successes within the sequence.<br>
!>  There are \f$\binom{n}{k}\f$ such sequences, since the binomial coefficient \f$\binom{n}{k}\f$ counts the
!>  number of ways to choose the positions of the \f$k\f$ successes among the \f$n\f$ trials.<br>
!>  The binomial distribution is concerned with the probability of obtaining any of these sequences,
!>  meaning the probability of obtaining one of them (\f$p^k q^{n-k}\f$) must be added \f$\binom{n}{k}\f$ times,
!>  hence,
!>  \f{equation}{
!>      \large
!>      \Pr(X=k) = \binom{n}{k} p^{k} (1-p)^{n-k} ~.
!>  \f}
!>
!>  Cumulative Distribution Function
!>  --------------------------------
!>
!>  The cumulative distribution function can be expressed as:
!>  \f{equation}{
!>      \large
!>      \up{CDF}(k; n, p) = \Pr(X\leq k) = \sum_{i=0}^{\lfloor k\rfloor}{n \choose i}p^{i}(1-p)^{n-i} ~,
!>  \f}
!>  where \f$\lfloor k\rfloor\f$ is the *floor* under \f$k\f$, i.e., the greatest integer less than or equal to \f$k\f$.<br>
!>
!>  It can also be represented in terms of the [regularized incomplete Beta function](@ref pm_mathBeta), as follows:<br>
!>  \f{eqnarray}{
!>      \large
!>      \up{CDF}(k;n,p)
!>      &=& \Pr(X\leq k) \\
!>      &=& I_{1-p}(n-k,k+1) \\
!>      &=& (n-k){n \choose k}\int _{0}^{1-p}t^{n-k-1}(1-t)^{k}\,dt ~.
!>  \f}
!>
!>  Sums of binomials
!>  -----------------
!>
!>  If \f$X \sim B(n, p)\f$ and \f$Y \sim B(m, p)\f$ are independent binomial variables with the same probability \f$p\f$,
!>  then \f$X + Y\f$ is again a binomial variable; its distribution is \f$Z = X + Y \sim B(n+m, p)\f$:<br>
!>  \f{eqnarray}{
!>      \large
!>      \up{P}(Z=k)
!>      &=& \sum_{i=0}^{k}\left[{\binom {n}{i}}p^{i}(1-p)^{n-i}\right]\left[{\binom {m}{k-i}}p^{k-i}(1-p)^{m-k+i}\right] \\
!>      &=& \binom{n+m}{k} p^{k}(1-p)^{n+m-k} ~.
!>  \f}
!>
!>  A Binomial distributed random variable \f$X \sim B(n, p)\f$ can be considered as the sum of \f$n\f$
!>  [Bernoulli distributed](@ref pm_distBern) random variables.<br>
!>  Thus, the sum of two Binomial distributed random variable \f$X \sim B(n, p)\f$ and \f$Y \sim B(m, p)\f$
!>  is equivalent to the sum of \f$n + m\f$ Bernoulli distributed random variables,
!>  which means \f$Z = X + Y \sim B(n + m, p)\f$.<br>
!>
!>  However, if \f$X\f$ and \f$Y\f$ do not have the same probability \f$p\f$,
!>  then the variance of the sum will be smaller than the variance of a binomial variable distributed as \f$B(n+m,{\bar{p}})\f$.<br>
!>
!>  \test
!>  [test_pm_distBinom](@ref test_pm_distBinom)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distBinom

    use pm_kind, only: SK, IK, LK, RKB
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distBinom"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Binomial
    !>  as defined in the description of [pm_distBinom](@ref pm_distBinom).
    !>
    !>  \details
    !>  See the documentation of [pm_distBinom](@ref pm_distBinom) for the definition of the Binomial distribution.
    !>
    !>  \interface{distBinom_type}
    !>  \code{.F90}
    !>
    !>      use pm_distBinom, only: distBinom_type
    !>      type(distBinom_type) :: distBinom
    !>
    !>      distBinom = distBinom_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distBinom](@ref test_pm_distBinom)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distBinom_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distBinom_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Binomial distribution</b> for an input `nsuc` within the discrete integer support of the distribution.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_distBinom](@ref pm_distBinom) for more information.<br>
    !>
    !>  \param[in]  nsuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values (numbers of success) at which the PMF must be computed.<br>
    !>  \param[in]  ntry        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the total number of trials in the Bernoulli Process.<br>
    !>  \param[in]  psuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the probability of success in the Bernoulli Process.<br>
    !>
    !>  \return
    !>  `logPMF`                :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `psuc`, containing the natural logarithm of the PMF
    !>                              of the distribution at the specified input `nsuc`.<br>
    !>
    !>  \interface{getBinomLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distBinom, only: getBinomLogPMF
    !>
    !>      logPMF = getBinomLogPMF(nsuc, ntry, psuc)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= psuc` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= nsuc` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= ntry` must hold for the corresponding input arguments.<br>
    !>  The condition `nsuc <= ntry` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setBinomLogPMF](@ref pm_distBinom::setBinomLogPMF)<br>
    !>
    !>  \example{getBinomLogPMF}
    !>  \include{lineno} example/pm_distBinom/getBinomLogPMF/main.F90
    !>  \compilef{getBinomLogPMF}
    !>  \output{getBinomLogPMF}
    !>  \include{lineno} example/pm_distBinom/getBinomLogPMF/main.out.F90
    !>  \postproc{getBinomLogPMF}
    !>  \include{lineno} example/pm_distBinom/getBinomLogPMF/main.py
    !>  \vis{getBinomLogPMF}
    !>  \image html pm_distBinom/getBinomLogPMF/getBinomLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBinom](@ref test_pm_distBinom)
    !>
    !>  \final{getBinomLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBinomLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBinomLogPMF_RK5(nsuc, ntry, psuc) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomLogPMF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBinomLogPMF_RK4(nsuc, ntry, psuc) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomLogPMF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBinomLogPMF_RK3(nsuc, ntry, psuc) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomLogPMF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBinomLogPMF_RK2(nsuc, ntry, psuc) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomLogPMF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: logPMF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBinomLogPMF_RK1(nsuc, ntry, psuc) result(logPMF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomLogPMF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: logPMF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Mass Function (PMF) of the
    !>  <b>Binomial distribution</b> for an input `nsuc` within the discrete integer support of the distribution.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_distBinom](@ref pm_distBinom) for more information.<br>
    !>
    !>  \param[out] logPMF      :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the natural logarithm of the PMF of the distribution at the specified `nsuc`.<br>
    !>  \param[in]  nsuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values (numbers of success) at which the PMF must be computed.<br>
    !>  \param[in]  ntry        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the total number of trials in the Bernoulli Process.<br>
    !>  \param[in]  logp        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as the output argument `logPMF`,
    !>                              containing the natural logarithm of the probability of success in the Bernoulli Process.<br>
    !>  \param[in]  logq        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as the output argument `logPMF`,
    !>                              containing the natural logarithm of the probability of failure in the Bernoulli Process.<br>
    !>
    !>  \interface{setBinomLogPMF}
    !>  \code{.F90}
    !>
    !>      use pm_distBinom, only: setBinomLogPMF
    !>
    !>      call setBinomLogPMF(logPMF, nsuc, ntry, logp, logq)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `logp <= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `logq <= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= nsuc` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= ntry` must hold for the corresponding input arguments.<br>
    !>  The condition `nsuc <= ntry` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(1 - exp(logp) - exp(logq)) <= epsilon(logp)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getBinomLogPMF](@ref pm_distBinom::getBinomLogPMF)<br>
    !>
    !>  \example{setBinomLogPMF}
    !>  \include{lineno} example/pm_distBinom/setBinomLogPMF/main.F90
    !>  \compilef{setBinomLogPMF}
    !>  \output{setBinomLogPMF}
    !>  \include{lineno} example/pm_distBinom/setBinomLogPMF/main.out.F90
    !>  \postproc{setBinomLogPMF}
    !>  \include{lineno} example/pm_distBinom/setBinomLogPMF/main.py
    !>  \vis{setBinomLogPMF}
    !>  \image html pm_distBinom/setBinomLogPMF/setBinomLogPMF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBinom](@ref test_pm_distBinom)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setBinomLogPMF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBinomLogPMF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBinomLogPMF_RK5(logPMF, nsuc, ntry, logp, logq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomLogPMF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPMF
        real(RKG)   , intent(in)                    :: logp, logq
        integer(IK) , intent(in)                    :: nsuc, ntry
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBinomLogPMF_RK4(logPMF, nsuc, ntry, logp, logq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomLogPMF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPMF
        real(RKG)   , intent(in)                    :: logp, logq
        integer(IK) , intent(in)                    :: nsuc, ntry
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBinomLogPMF_RK3(logPMF, nsuc, ntry, logp, logq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomLogPMF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPMF
        real(RKG)   , intent(in)                    :: logp, logq
        integer(IK) , intent(in)                    :: nsuc, ntry
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBinomLogPMF_RK2(logPMF, nsuc, ntry, logp, logq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomLogPMF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPMF
        real(RKG)   , intent(in)                    :: logp, logq
        integer(IK) , intent(in)                    :: nsuc, ntry
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBinomLogPMF_RK1(logPMF, nsuc, ntry, logp, logq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomLogPMF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPMF
        real(RKG)   , intent(in)                    :: logp, logq
        integer(IK) , intent(in)                    :: nsuc, ntry
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  <b>Binomial distribution</b> for an input `nsuc` within the discrete integer support of the distribution.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_distBinom](@ref pm_distBinom) for more information.<br>
    !>
    !>  \param[in]  nsuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values (numbers of success) at which the CDF must be computed.<br>
    !>  \param[in]  ntry        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the total number of trials in the Bernoulli Process.<br>
    !>  \param[in]  psuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the probability of success in the Bernoulli Process.<br>
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar (or array of the same rank, shape, and size as other array-like arguments),
    !>                              of the same type and kind as `logp`, containing the CDF of
    !>                              the distribution at the specified point.<br>
    !>
    !>  \interface{getBinomCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distBinom, only: getBinomCDF
    !>
    !>      cdf = getBinomCDF(nsuc, ntry, psuc)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= nsuc` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= ntry` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= psuc` must hold for the corresponding input arguments.<br>
    !>  The condition `psuc <= 1` must hold for the corresponding input arguments.<br>
    !>  The condition `nsuc <= ntry` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setBinomCDF](@ref pm_distBinom::setBinomCDF)<br>
    !>
    !>  \example{getBinomCDF}
    !>  \include{lineno} example/pm_distBinom/getBinomCDF/main.F90
    !>  \compilef{getBinomCDF}
    !>  \output{getBinomCDF}
    !>  \include{lineno} example/pm_distBinom/getBinomCDF/main.out.F90
    !>  \postproc{getBinomCDF}
    !>  \include{lineno} example/pm_distBinom/getBinomCDF/main.py
    !>  \vis{getBinomCDF}
    !>  \image html pm_distBinom/getBinomCDF/getBinomCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBinom](@ref test_pm_distBinom)
    !>
    !>  \final{getBinomCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBinomCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBinomCDF_RK5(nsuc, ntry, psuc) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBinomCDF_RK4(nsuc, ntry, psuc) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBinomCDF_RK3(nsuc, ntry, psuc) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBinomCDF_RK2(nsuc, ntry, psuc) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBinomCDF_RK1(nsuc, ntry, psuc) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinomCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the <b>Binomial distribution</b>.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_distBinom](@ref pm_distBinom) for more information.<br>
    !>
    !>  \param[out] cdf         :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the distribution CDF at the specified `nsuc`.<br>
    !>  \param[in]  nsuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the values (numbers of success) at which the CDF must be computed.<br>
    !>  \param[in]  ntry        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              type `integer` of default kind \IK, containing the total number of trials in the Bernoulli Process.<br>
    !>  \param[in]  psuc        :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the probability of success in the Bernoulli Process.<br>
    !>  \param[out] info        :   The output scalar of type `integer` of default kind \IK.<br>
    !>                              On output, it is zero if and only if the CDF computation converges.<br>
    !>                              For more information, see the documentation of the corresponding argument of [setBetaInc](@ref pm_mathBeta::setBetaInc).<br>
    !>
    !>  \interface{setBinomCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distBinom, only: setBinomCDF
    !>
    !>      call setBinomCDF(cdf, nsuc, ntry, psuc, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= nsuc` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= ntry` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= psuc` must hold for the corresponding input arguments.<br>
    !>  The condition `psuc <= 1` must hold for the corresponding input arguments.<br>
    !>  The condition `nsuc <= ntry` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getBinomCDF](@ref pm_distBinom::getBinomCDF)<br>
    !>
    !>  \example{setBinomCDF}
    !>  \include{lineno} example/pm_distBinom/setBinomCDF/main.F90
    !>  \compilef{setBinomCDF}
    !>  \output{setBinomCDF}
    !>  \include{lineno} example/pm_distBinom/setBinomCDF/main.out.F90
    !>  \postproc{setBinomCDF}
    !>  \include{lineno} example/pm_distBinom/setBinomCDF/main.py
    !>  \vis{setBinomCDF}
    !>  \image html pm_distBinom/setBinomCDF/setBinomCDF.IK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBinom](@ref test_pm_distBinom)
    !>
    !>  \final{setBinomCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBinomCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBinomCDF_RK5(cdf, nsuc, ntry, psuc, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBinomCDF_RK4(cdf, nsuc, ntry, psuc, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBinomCDF_RK3(cdf, nsuc, ntry, psuc, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBinomCDF_RK2(cdf, nsuc, ntry, psuc, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBinomCDF_RK1(cdf, nsuc, ntry, psuc, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBinomCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: nsuc, ntry
        real(RKG)   , intent(in)                    :: psuc
        real(RKG)   , intent(out)                   :: cdf
        integer(IK) , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distBinom
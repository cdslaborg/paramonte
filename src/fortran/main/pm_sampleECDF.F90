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
!>  This module contains classes and procedures for computing the Empirical Cumulative Distribution Function (ECDF)
!>  of an observational sample and the associated the various properties.
!>
!>  \details
!>  An **empirical Cumulative Distribution Function** (eCDF) is the distribution function associated with the empirical measure of a sample.<br>
!>  This cumulative distribution function is a step function that jumps up by \f$1 / N\f$ at each of the \f$N\f$ data points.<br>
!>  Its value at any specified value of the measured variable is the fraction of observations of the measured variable that are less than or equal to the specified value.<br>
!>  The empirical distribution function is an estimate of the cumulative distribution function that generated the points in the sample.<br>
!>  It converges with probability 1 to that underlying distribution, according to the Glivenko–Cantelli theorem.<br>
!>  A number of results exist to quantify the rate of convergence of the empirical distribution function to the underlying cumulative distribution function.<br>
!>
!>  **Definition**<br>
!>
!>  Let \f$(X_1, \ldots, X_N)\f$ be independent, identically distributed real random variables with the common cumulative distribution function \f$F(t)\f$.<br>
!>  Then the empirical distribution function is defined as,
!>  \f{equation}{
!>      {\widehat{F}}_{N}(t) = \frac{{\mbox{number of elements in the sample}}\leq t}{N} = {\frac{1}{N}} \sum_{i = 1}^{N} \mathbf{1}_{X_{i}\leq t} ~,
!>  \f}
!>  where \f${\mathbf{1}}_{{A}}\f$ is the indicator of event \f$A\f$.<br>
!>  For a fixed \f$t\f$, the indicator \f$\mathbf{1}_{X_{i}\leq t}\f$ is a Bernoulli random variable with parameter \f$p = F(t)\f$.<br>
!>  Hence, \f$N{\widehat{F}}_{N}(t)\f$ is a **binomial random variable** with mean \f$N\times F(t)\f$ and variance \f$N\times F(t)(1 − F(t))\f$.<br>
!>  This implies that \f${\widehat{F}}_{N}(t)\f$ is an unbiased estimator for \f$F(t)\f$.<br>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>  [pm_sampleACT](@ref pm_sampleACT)<br>
!>  [pm_sampleCCF](@ref pm_sampleCCF)<br>
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>  [pm_sampleCov](@ref pm_sampleCov)<br>
!>  [pm_sampleConv](@ref pm_sampleConv)<br>
!>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
!>  [pm_sampleMean](@ref pm_sampleMean)<br>
!>  [pm_sampleNorm](@ref pm_sampleNorm)<br>
!>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
!>  [pm_sampleScale](@ref pm_sampleScale)<br>
!>  [pm_sampleShift](@ref pm_sampleShift)<br>
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_sampleAffinity](@ref pm_sampleAffinity)<br>
!>  [pm_sampleVar](@ref pm_sampleVar)<br>
!>  [Empirical distribution function](https://en.wikipedia.org/wiki/Empirical_distribution_function)<br>
!>
!>  \test
!>  [test_pm_sampleECDF](@ref test_pm_sampleECDF)<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleECDF

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleECDF"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return the Empirical Cumulative Distribution Function (ECDF) of a univariate (optionally weighted) sample of size `size(ecdf)`.
    !>
    !>  \param[out]     ecdf    :   The output `contiguous` array of rank `1` of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              of the same size `nsam` as the number of observations in the target sample whose ECDF must be computed.<br>
    !>                              On output, it contains the Empirical Cumulative Distribution Function (ECDF) of the input sample.<br>
    !>  \param[in]      weight  :   The input `contiguous` vector of length `nsam` of,
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK, or
    !>                                  <li>    type `real` of the same kind as the kind of `ecdf`,
    !>                              </ol>
    !>                              containing the corresponding weights of individual `nsam` observations in the target sample.<br>
    !>                              (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled).)
    !>  \param[in]      weisum  :   The input scalar of the same type and kind as `weight` containing the quantity `sum(weight)`.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `weight` is also present.)
    !>  \param[out]     lcdf    :   The output `contiguous` array of the same type, kind, rank, and size as the output `ecdf`
    !>                              containing the lower confidence bound on the ECDF as the specified \f$1 - \alpha\f$ probability.<br>
    !>                              (**optional**, if missing, it is not computed. It can be present **if and only if** the input argument `weight` is missing or is of type `integer` of default kind \IK.)
    !>  \param[out]     ucdf    :   The output `contiguous` array of the same type, kind, rank, and size as the output `ecdf`
    !>                              containing the upper confidence bound on the ECDF as the specified \f$1 - \alpha\f$ probability.<br>
    !>                              (**optional**, if missing, it is not computed. It can be present **if and only if** the input argument `weight` is missing or is of type `integer` of default kind \IK.)
    !>  \param[in]      alpha   :   The input scalar of type `real` of the same kind as the kind of `ecdf` such that \f$1 - \alpha\f$ represents
    !>                              the probability of the parent CDF being bounded by the output upper and lower confidence bounds of the output ecdf.<br>
    !>                              For example, if \f$\alpha = 0.05\f$, then the output `lcdf` and `ucdf` confidence bounds contain
    !>                              the CDF of the true parent distribution of the sample with \f$ 95\% \f$ probability.<br>
    !>                              In other words,<br>
    !>                              <ol>
    !>                                  <li>    the lower bound marks the \f$\alpha / 2\f$ probability contour.<br>
    !>                                  <li>    the upper bound marks the \f$1 - \alpha / 2\f$ probability contour.<br>
    !>                              </ol>
    !>                              (**optional**, default = `0.05`. It can be present **if and only if** the input argument `weight` is missing or is of type `integer` of default kind \IK.)
    !>
    !>  \interface{setECDF}
    !>  \code{.F90}
    !>
    !>      use pm_sampleECDF, only: setECDF
    !>
    !>      call setECDF(ecdf(1 : nsam)                  , lcdf = lcdf(1 : nsam), ucdf = ucdf(1 : nsam), alpha = alpha)
    !>      call setECDF(ecdf(1 : nsam), weight(1 : nsam), lcdf = lcdf(1 : nsam), ucdf = ucdf(1 : nsam), alpha = alpha) ! integer (frequency) weight
    !>      call setECDF(ecdf(1 : nsam), weight(1 : nsam)) ! real (reliability) weight
    !>      !
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [setCumSum](@ref pm_mathCumSum::setCumSum)<br>
    !>  [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp)<br>
    !>  [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_sampleECDF/setECDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_sampleECDF/setECDF/main.out.F90
    !>  \vis
    !>  \image html pm_sampleECDF/setECDF/main.norm.10.out.png width=700
    !>  \image html pm_sampleECDF/setECDF/main.norm.100.out.png width=700
    !>  \image html pm_sampleECDF/setECDF/main.norm.1000.out.png width=700
    !>  \image html pm_sampleECDF/setECDF/main.norm.10000.out.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleECDF](@ref test_pm_sampleECDF)
    !>
    !>  \final{setECDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface setECDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setECDF_ONE_D1_RK5(ecdf, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_ONE_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setECDF_ONE_D1_RK4(ecdf, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_ONE_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setECDF_ONE_D1_RK3(ecdf, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_ONE_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setECDF_ONE_D1_RK2(ecdf, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_ONE_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setECDF_ONE_D1_RK1(ecdf, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_ONE_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setECDF_WIK_D1_RK5(ecdf, weight, weisum, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WIK_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK) , intent(in)                            :: weisum
        integer(IK) , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setECDF_WIK_D1_RK4(ecdf, weight, weisum, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WIK_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK) , intent(in)                            :: weisum
        integer(IK) , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setECDF_WIK_D1_RK3(ecdf, weight, weisum, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WIK_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK) , intent(in)                            :: weisum
        integer(IK) , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setECDF_WIK_D1_RK2(ecdf, weight, weisum, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WIK_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK) , intent(in)                            :: weisum
        integer(IK) , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setECDF_WIK_D1_RK1(ecdf, weight, weisum, lcdf, ucdf, alpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WIK_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK) , intent(in)                            :: weisum
        integer(IK) , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: ucdf(:)
        real(TKG)   , intent(out)   , contiguous, optional  :: lcdf(:)
        real(TKG)   , intent(in)                , optional  :: alpha
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setECDF_WRK_D1_RK5(ecdf, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WRK_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)   , intent(in)                            :: weisum
        real(TKG)   , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setECDF_WRK_D1_RK4(ecdf, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WRK_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)   , intent(in)                            :: weisum
        real(TKG)   , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setECDF_WRK_D1_RK3(ecdf, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WRK_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)   , intent(in)                            :: weisum
        real(TKG)   , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setECDF_WRK_D1_RK2(ecdf, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WRK_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)   , intent(in)                            :: weisum
        real(TKG)   , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setECDF_WRK_D1_RK1(ecdf, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setECDF_WRK_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)   , intent(in)                            :: weisum
        real(TKG)   , intent(in)    , contiguous            :: weight(:)
        real(TKG)   , intent(out)   , contiguous            :: ecdf(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleECDF ! LCOV_EXCL_LINE
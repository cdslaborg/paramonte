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
!>  This module contains classes and procedures for performing various statistical tests.
!>
!>  \details
!>
!>  Kolmogorov-Smirnov (KS) Test
!>  ----------------------------
!>
!>  The Kolmogorov–Smirnov test (K–S test or KS test) is a **nonparametric test** of the equality of continuous, one-dimensional probability distributions
!>  that can be used to compare a sample with a reference probability distribution (one-sample K–S test), or to compare two samples (two-sample K–S test).<br>
!>  In essence, the test answers the question *How likely is it that we would see a collection of samples like this if they were drawn from that probability distribution?* or,
!>  in the second case, *How likely is it that we would see two sets of samples like this if they were drawn from the same (but unknown) probability distribution?*.<br>
!>  It is named after **Andrey Kolmogorov** and **Nikolai Smirnov**.<br>
!>
!>  The Kolmogorov–Smirnov statistic quantifies a [distance](@ref pm_distanceKolm) between the [empirical distribution function](@ref pm_sampleECDF) of
!>  the sample and the cumulative distribution function of the reference distribution, or between the empirical distribution functions of two samples.<br>
!>  The **null distribution** of this statistic is calculated under the **null hypothesis** that the **sample is drawn from the reference distribution**
!>  (in the one-sample case) or that the samples are drawn from the same distribution (in the two-sample case).<br>
!>  In the **one-sample** case, the distribution considered under the null hypothesis may be continuous, purely discrete or mixed.<br>
!>  In the **two-sample** case, the distribution considered under the null hypothesis is a continuous distribution but is otherwise unrestricted.<br>
!>  However, the two sample test can also be performed under more general conditions that allow for discontinuity, heterogeneity and dependence across samples.<br>
!>  The two-sample K–S test is one of the most useful and general nonparametric methods for **comparing two samples**, as it is **sensitive**
!>  to differences in both **location** and **shape** of the empirical cumulative distribution functions of the two samples.<br>
!>  The Kolmogorov–Smirnov test can be modified to serve as a **goodness of fit test**.<br>
!>
!>  \note
!>  The **effective sample size** in a two-sample KS test is computed as,
!>  \f{equation}
!>      N_e = \frac{N_1 N_2}{N_1 + N_2} ~,
!>  \f}
!>  where \f$N_1\f$ and \f$N_2\f$ represent the sizes of the first and the second samples in the test respectively.<br>
!>  The nature of the approximation involved in the two-sample KS test is that it becomes asymptotically accurate as
!>  the effective sample size \f$N_e\f$ becomes large, but the approximation is reasonable even for as low as \f$N_e\approx 4\f$.<br>
!>
!>  Testing for Normality
!>  =====================
!>
!>  In the special case of testing for normality of the distribution, samples are [shifted](@ref pm_sampleShift),
!>  [scaled](@ref pm_sampleScale), or [standardized](@ref pm_sampleNorm) and compared with a standard normal distribution.<br>
!>  This is equivalent to setting the mean and variance of the reference distribution equal to the sample estimates,
!>  and it is known that using these to define the specific reference distribution changes the null distribution of the test statistic.<br>
!>  Various studies have found that, even in this corrected form, the test is less powerful for testing normality than the **Shapiro–Wilk test** or **Anderson–Darling test**.<br>
!>  However, these other tests have their own disadvantages.<br>
!>  For instance the Shapiro–Wilk test is known not to work well in samples with many identical values.<br>
!>
!>  \see
!>  [pm_statest](@ref pm_statest)<br>
!>  [pm_distKolm](@ref pm_distKolm)<br>
!>  [pm_distanceKolm](@ref pm_distanceKolm)<br>
!>
!>  \test
!>  [test_pm_statest](@ref test_pm_statest)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_statest

    use pm_kind, only: IK, RK, SK
    use pm_distanceKolm, only: getDisKolm, setDisKolm
    use pm_distanceKolm, only: ascending, ascending_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_statest"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the probability of the null-hypothesis that `sample1` of size `nsam1` originates from the same
    !>  distribution as that of `sample2` of size `nsam2` or from the Uniform distribution or other distribution whose custom CDF is given.<br>
    !>
    !>  \details
    !>  See [pm_statest](@ref pm_statest) for the mathematical definition of the KS test.<br>
    !>
    !>  \param[in]      statKS      :   The input scalar of the same type and kind as the output `probKS`,
    !>                                  representing the KS test statistic for the null-hypothesis considered.<br>
    !>                                  This quantity is the same as the [Kolmogorov distance](@ref pm_distanceKolm)
    !>                                  and can be readily obtained for any two samples or a sample against a distribution CDF
    !>                                  via the generic interfaces [getDisKolm](@ref pm_distanceKolm::getDisKolm) or
    !>                                  [setDisKolm](@ref pm_distanceKolm::setDisKolm).<br>
    !>  \param[in]      weisum1     :   The input scalar of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, <br>
    !>                                      <li>    type `real` of the same kind as that of the output `probKS` of type `real`, <br>
    !>                                  </ol>
    !>                                  representing either,
    !>                                  <ol>
    !>                                      <li>    the size of the first **unweighted** sample (if `weisum1` is of type `integer`) or,
    !>                                      <li>    the quantity `sum(weight1)` where `weight1` is the vector of weights of the first sample in the KS test,
    !>                                  </ol>
    !>  \param[in]      weisum2     :   The input scalar of the same type and kind as the input `weisum1`,
    !>                                  representing either,
    !>                                  <ol>
    !>                                      <li>    the size of the second **unweighted** sample (if `weisum2` is of type `integer`) or,
    !>                                      <li>    the quantity `sum(weight2)` where `weight2` is the vector of weights of the second sample in the KS test,
    !>                                  </ol>
    !>                                  (**optional**. It must be present if and only if the input argument `weisum1` is also present and the KS test involves two samples.)
    !>  \param[in]      wsqsum1     :   The input scalar of type `real` of the same kind as that of the input `weisum1` of type `real`,
    !>                                  representing the quantity `sum(weight1**2)` where `weight1` is the vector of weights of the first sample in the KS test.<br>
    !>                                  This quantity must be supplied **if and only if** the sample weights are [reliability weights](@ref pm_sampleWeight),
    !>                                  which requires the weights (and hence, `weisum1` and `wsqsum1`) to be of type `real`.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `weisum1` is present and is of type `real`.)
    !>  \param[in]      wsqsum2     :   The input scalar of type `real` of the same kind as that of the input `weisum2` of type `real`,
    !>                                  representing the quantity `sum(weight2**2)` where `weight2` is the vector of weights of the first sample in the KS test.<br>
    !>                                  This quantity must be supplied **if and only if** the sample weights are [reliability weights](@ref pm_sampleWeight),
    !>                                  which requires the weights (and hence, `weisum2` and `wsqsum2`) to be of type `real`.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `weisum2` is present and is of type `real`.)
    !>
    !>  \return
    !>  `probKS`                    :   The output scalar of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  representing the probability of observing a KS test statistic as extreme
    !>                                  as or more extreme than the observed value under the null hypothesis.<br>
    !>                                  Small values of `probKS` cast doubt on the validity of the null hypothesis.<br>
    !>                                  In other words, `probKS` represents the probability that the specified two samples
    !>                                  or the sample and the specified CDF originate from the same distribution.<br>
    !>
    !>  \interface{getProbKS}
    !>  \code{.F90}
    !>
    !>      use pm_statest, only: getProbKS, ascending
    !>
    !>      ! one-sample KS test.
    !>
    !>      probKS = getProbKS(statKS, weisum1) ! only unweighted or (integer) frequency-weighted sample.
    !>      probKS = getProbKS(statKS, weisum1, wsqsum1) ! only (real) reliability-weighted sample.
    !>
    !>      ! two-sample KS test.
    !>
    !>      probKS = getProbKS(statKS, weisum1, weisum2) ! only unweighted or (integer) frequency-weighted samples.
    !>      probKS = getProbKS(statKS, weisum1, weisum2, wsqsum1) ! only if `sample1` is reliability weighted and `sample2` is unweighted or frequency-weighted.
    !>      probKS = getProbKS(statKS, weisum1, weisum2, wsqsum1, wsqsum2) ! only if both samples are reliability-real-weighted samples.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < weisum1` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < weisum2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < wsqsum1` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < wsqsum2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= statKS` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \example{getProbKS}
    !>  \include{lineno} example/pm_statest/getProbKS/main.F90
    !>  \compilef{getProbKS}
    !>  \output{getProbKS}
    !>  \include{lineno} example/pm_statest/getProbKS/main.out.F90
    !>
    !>  \test
    !>  [test_pm_statest](@ref test_pm_statest)<br>
    !>
    !>  \naming
    !>  \code{.F90}
    !>      getProbKS_WIX_D0_RK5()
    !>                ||| || |||
    !>                ||| || |||
    !>                ||| || |||
    !>                ||| || |The Kind of the output.
    !>                ||| The rank of the input arguments.
    !>                The sample presence and weight types: X => missing, WDD => default (unweighted) / default (unweighted), WID => integer-weighted, default, WRD => real-weighted, default., WII => integer-weighted, integer-weighted, WRR => real-weighted, real-weighted.
    !>  \endcode
    !>
    !>  \final{getProbKS}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! SX

    interface getProbKS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getProbKS_WIX_D0_RK5(statKS, weisum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WIX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK4_ENABLED
    PURE module function getProbKS_WIX_D0_RK4(statKS, weisum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WIX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK3_ENABLED
    PURE module function getProbKS_WIX_D0_RK3(statKS, weisum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WIX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK2_ENABLED
    PURE module function getProbKS_WIX_D0_RK2(statKS, weisum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WIX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK1_ENABLED
    PURE module function getProbKS_WIX_D0_RK1(statKS, weisum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WIX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getProbKS_WRX_D0_RK5(statKS, weisum1, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK4_ENABLED
    PURE module function getProbKS_WRX_D0_RK4(statKS, weisum1, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK3_ENABLED
    PURE module function getProbKS_WRX_D0_RK3(statKS, weisum1, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK2_ENABLED
    PURE module function getProbKS_WRX_D0_RK2(statKS, weisum1, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK1_ENABLED
    PURE module function getProbKS_WRX_D0_RK1(statKS, weisum1, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS

    interface getProbKS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getProbKS_WII_D0_RK5(statKS, weisum1, weisum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WII_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK4_ENABLED
    PURE module function getProbKS_WII_D0_RK4(statKS, weisum1, weisum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WII_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK3_ENABLED
    PURE module function getProbKS_WII_D0_RK3(statKS, weisum1, weisum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WII_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK2_ENABLED
    PURE module function getProbKS_WII_D0_RK2(statKS, weisum1, weisum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WII_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK1_ENABLED
    PURE module function getProbKS_WII_D0_RK1(statKS, weisum1, weisum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WII_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getProbKS_WRI_D0_RK5(statKS, weisum1, weisum2, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRI_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK4_ENABLED
    PURE module function getProbKS_WRI_D0_RK4(statKS, weisum1, weisum2, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRI_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK3_ENABLED
    PURE module function getProbKS_WRI_D0_RK3(statKS, weisum1, weisum2, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRI_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK2_ENABLED
    PURE module function getProbKS_WRI_D0_RK2(statKS, weisum1, weisum2, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRI_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK1_ENABLED
    PURE module function getProbKS_WRI_D0_RK1(statKS, weisum1, weisum2, wsqsum1) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRI_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getProbKS_WRR_D0_RK5(statKS, weisum1, weisum2, wsqsum1, wsqsum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK4_ENABLED
    PURE module function getProbKS_WRR_D0_RK4(statKS, weisum1, weisum2, wsqsum1, wsqsum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK3_ENABLED
    PURE module function getProbKS_WRR_D0_RK3(statKS, weisum1, weisum2, wsqsum1, wsqsum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK2_ENABLED
    PURE module function getProbKS_WRR_D0_RK2(statKS, weisum1, weisum2, wsqsum1, wsqsum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

#if RK1_ENABLED
    PURE module function getProbKS_WRR_D0_RK1(statKS, weisum1, weisum2, wsqsum1, wsqsum2) result(probKS)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS_WRR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)                               :: probKS
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the probability and the corresponding [Kolmogorov distribution](@ref pm_distKolm) quantile of the null-hypothesis that `sample1`
    !>  of size `nsam1` originates from the same distribution as that of `sample2` of size `nsam2` or
    !>  from the Uniform distribution or other distribution whose custom CDF is given.<br>
    !>
    !>  \details
    !>  See [pm_statest](@ref pm_statest) for the mathematical definition of the KS test.<br>
    !>
    !>  \param[out]     probKS      :   The output scalar of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  representing the probability of observing a KS test statistic as extreme
    !>                                  as or more extreme than the observed value under the null hypothesis.<br>
    !>                                  Small values of `probKS` cast doubt on the validity of the null hypothesis.<br>
    !>                                  In other words, `probKS` represents the probability that the specified two samples
    !>                                  or the sample and the specified CDF originate from the same distribution.<br>
    !>  \param[out]     quanKS      :   The output scalar of the same type and kind as the output `probKS`,
    !>                                  containing the [Kolmogorov distribution](@ref pm_distKolm) quantile
    !>                                  corresponding to the input KS statistic `statKS` and sample size(s).<br>
    !>  \param[in]      statKS      :   The input scalar of the same type and kind as the output `probKS`,
    !>                                  representing the KS test statistic for the null-hypothesis considered.<br>
    !>                                  This quantity is the same as the [Kolmogorov distance](@ref pm_distanceKolm)
    !>                                  and can be readily obtained for any two samples or a sample against a distribution CDF
    !>                                  via the generic interfaces [getDisKolm](@ref pm_distanceKolm::getDisKolm) or
    !>                                  [setDisKolm](@ref pm_distanceKolm::setDisKolm).<br>
    !>  \param[in]      weisum1     :   The input scalar of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, <br>
    !>                                      <li>    type `real` of the same kind as that of the output `probKS` of type `real`, <br>
    !>                                  </ol>
    !>                                  representing either,
    !>                                  <ol>
    !>                                      <li>    the size of the first **unweighted** sample (if `weisum1` is of type `integer`) or,
    !>                                      <li>    the quantity `sum(weight1)` where `weight1` is the vector of weights of the first sample in the KS test,
    !>                                  </ol>
    !>  \param[in]      weisum2     :   The input scalar of the same type and kind as the input `weisum1`,
    !>                                  representing either,
    !>                                  <ol>
    !>                                      <li>    the size of the second **unweighted** sample (if `weisum2` is of type `integer`) or,
    !>                                      <li>    the quantity `sum(weight2)` where `weight2` is the vector of weights of the second sample in the KS test,
    !>                                  </ol>
    !>                                  (**optional**. It must be present if and only if the input argument `weisum1` is also present and the KS test involves two samples.)
    !>  \param[in]      wsqsum1     :   The input scalar of type `real` of the same kind as that of the input `weisum1` of type `real`,
    !>                                  representing the quantity `sum(weight1**2)` where `weight1` is the vector of weights of the first sample in the KS test.<br>
    !>                                  This quantity must be supplied **if and only if** the sample weights are [reliability weights](@ref pm_sampleWeight),
    !>                                  which requires the weights (and hence, `weisum1` and `wsqsum1`) to be of type `real`.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `weisum1` is present and is of type `real`.)
    !>  \param[in]      wsqsum2     :   The input scalar of type `real` of the same kind as that of the input `weisum2` of type `real`,
    !>                                  representing the quantity `sum(weight2**2)` where `weight2` is the vector of weights of the first sample in the KS test.<br>
    !>                                  This quantity must be supplied **if and only if** the sample weights are [reliability weights](@ref pm_sampleWeight),
    !>                                  which requires the weights (and hence, `weisum2` and `wsqsum2`) to be of type `real`.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `weisum2` is present and is of type `real`.)
    !>
    !>  \interface{setProbKS}
    !>  \code{.F90}
    !>
    !>      use pm_statest, only: setProbKS, ascending
    !>
    !>      ! one-sample KS test.
    !>
    !>      call setProbKS(probKS, quanKS, statKS, weisum1) ! only unweighted or (integer) frequency-weighted sample.
    !>      call setProbKS(probKS, quanKS, statKS, weisum1, wsqsum1) ! only (real) reliability-weighted sample.
    !>
    !>      ! two-sample KS test.
    !>
    !>      call setProbKS(probKS, quanKS, statKS, weisum1, weisum2) ! only unweighted or (integer) frequency-weighted samples.
    !>      call setProbKS(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1) ! only if `sample1` is reliability weighted and `sample2` is unweighted or frequency-weighted.
    !>      call setProbKS(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1, wsqsum2) ! only if both samples are reliability-real-weighted samples.
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < weisum1` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < weisum2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < wsqsum1` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < wsqsum2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= statKS` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \example{setProbKS}
    !>  \include{lineno} example/pm_statest/setProbKS/main.F90
    !>  \compilef{setProbKS}
    !>  \output{setProbKS}
    !>  \include{lineno} example/pm_statest/setProbKS/main.out.F90
    !>
    !>  \test
    !>  [test_pm_statest](@ref test_pm_statest)<br>
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setProbKS_WIX_D0_RK5()
    !>                ||| || |||
    !>                ||| || |||
    !>                ||| || |||
    !>                ||| || |The Kind of the output.
    !>                ||| The rank of the input arguments.
    !>                The sample presence and weight types: X => missing, WDD => default (unweighted) / default (unweighted), WID => integer-weighted, default, WRD => real-weighted, default., WII => integer-weighted, integer-weighted, WRR => real-weighted, real-weighted.
    !>  \endcode
    !>
    !>  \final{setProbKS}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! SX

    interface setProbKS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setProbKS_WIX_D0_RK5(probKS, quanKS, statKS, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WIX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setProbKS_WIX_D0_RK4(probKS, quanKS, statKS, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WIX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setProbKS_WIX_D0_RK3(probKS, quanKS, statKS, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WIX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setProbKS_WIX_D0_RK2(probKS, quanKS, statKS, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WIX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setProbKS_WIX_D0_RK1(probKS, quanKS, statKS, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WIX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                :: weisum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setProbKS_WRX_D0_RK5(probKS, quanKS, statKS, weisum1, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRX_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setProbKS_WRX_D0_RK4(probKS, quanKS, statKS, weisum1, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRX_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setProbKS_WRX_D0_RK3(probKS, quanKS, statKS, weisum1, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRX_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setProbKS_WRX_D0_RK2(probKS, quanKS, statKS, weisum1, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRX_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setProbKS_WRX_D0_RK1(probKS, quanKS, statKS, weisum1, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRX_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS

    interface setProbKS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setProbKS_WII_D0_RK5(probKS, quanKS, statKS, weisum1, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WII_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setProbKS_WII_D0_RK4(probKS, quanKS, statKS, weisum1, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WII_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setProbKS_WII_D0_RK3(probKS, quanKS, statKS, weisum1, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WII_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setProbKS_WII_D0_RK2(probKS, quanKS, statKS, weisum1, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WII_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setProbKS_WII_D0_RK1(probKS, quanKS, statKS, weisum1, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WII_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                :: weisum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setProbKS_WRI_D0_RK5(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRI_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setProbKS_WRI_D0_RK4(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRI_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setProbKS_WRI_D0_RK3(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRI_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setProbKS_WRI_D0_RK2(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRI_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setProbKS_WRI_D0_RK1(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRI_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        integer(IK) , intent(in)                :: weisum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setProbKS_WRR_D0_RK5(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1, wsqsum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setProbKS_WRR_D0_RK4(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1, wsqsum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setProbKS_WRR_D0_RK3(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1, wsqsum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setProbKS_WRR_D0_RK2(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1, wsqsum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setProbKS_WRR_D0_RK1(probKS, quanKS, statKS, weisum1, weisum2, wsqsum1, wsqsum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProbKS_WRR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: weisum1, wsqsum1
        real(RKG)   , intent(in)                :: weisum2, wsqsum2
        real(RKG)   , intent(in)                :: statKS
        real(RKG)   , intent(out)               :: probKS
        real(RKG)   , intent(out)               :: quanKS
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_statest ! LCOV_EXCL_LINE
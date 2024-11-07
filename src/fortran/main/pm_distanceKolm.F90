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
!>  This module contains classes and procedures for computing the Kolmogorov statistical distance.
!>
!>  \details
!>  The Kolmogorov distance of a univariate observational sample from another univariate observational sample
!>  is the largest separation between the [Empirical Distribution Functions](@ref pm_sampleECDF) of the two samples.<br>
!>  Formally, the empirical distribution function \f$F_n\f$ for \f$n\f$ independent and identically distributed (i.i.d.)
!>  ordered observations \f$X_i\f$ is defined as,<br>
!>  \f{equation}{
!>      F_{n}(x) = {\frac{{\text{number of (elements in the sample}} \leq x)}{n}} = {\frac{1}{n}} \sum_{i=1}^{n}1_{(-\infty ,x]}(X_{i}) ~,
!>  \f}
!>  where \f$1_{(-\infty ,x]}(X_{i})\f$ is the indicator function, equal to \f$1\f$ if \f$X_{i}\leq x\f$ and equal to \f$0\f$ otherwise.<br>
!>  The Kolmogorov–Smirnov distance (or statistic) for a given cumulative distribution function \f$F(x)\f$ is,
!>  \f{equation}{
!>      D_{n} = \sup_{x}|F_{n}(x) - F(x)| ~,
!>  \f}
!>  where \f$\sup_x\f$ is the supremum of the set of distances.<br>
!>  Intuitively, the statistic takes the largest absolute difference
!>  between the two distribution functions across all \f$x\f$ values.<br>
!>  By the Glivenko–Cantelli theorem, if the sample comes from distribution \f$F(x)\f$,
!>  then \f$D_n\f$ converges to \f$0\f$ almost surely in the limit when \f$n\f$ goes to infinity.<br>
!>  Kolmogorov strengthened this result, by effectively providing the rate of this convergence
!>  through the definition of the [Kolmogorov distribution](@ref pm_distKolm).
!>
!>  \see
!>  [pm_distKolm](@ref pm_distKolm)<br>
!>
!>  \test
!>  [test_pm_distanceKolm](@ref test_pm_distanceKolm)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distanceKolm

    use pm_kind, only: IK, RK, SK
    use pm_array, only: ascending, ascending_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distanceKolm"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Kolmogorov distance of a `sample1` of size `nsam1` from another
    !>  sample `sample2` of size `nsam2` or the CDF of the Uniform or a custom reference distribution.<br>
    !>
    !>  \details
    !>  See [pm_distanceKolm](@ref pm_distanceKolm) for the mathematical definition of the Kolmogorov distance.
    !>
    !>  \param[in]      sample1     :   The input vector of size `nsam1` of the same type and kind as the output `disKolm`,
    !>                                  representing the first sample whose distance form the other sample must be computed.<br>
    !>  \param[in]      weight1     :   The input vector of the same size as the size of `sample1` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, <br>
    !>                                      <li>    type `real` of the same kind as that of the output `disKolm` of type `real`, <br>
    !>                                  </ol>
    !>                                  representing the weights of the corresponding elements of the input `sample1`.<br>
    !>                                  (**optional**. default = `[(1, i = 1, size(sample1))]`. It must be present if the input argument `weight2` is also present.)
    !>  \param[in]      weisum1     :   The input scalar of the same type and kind as the input `weight1`, containing the value `sum(weight1)`.<br>
    !>                                  (**optional**. It must be present if and only if the input argument `weight1` is also present.)
    !>  \param[in]      sample2     :   The input vector of size `nsam2` of the same type and kind as the input `disKolm`,
    !>                                  representing the second sample whose distance form the other sample must be computed.<br>
    !>                                  (**optional**, the default is the [Uniform distribution](@ref pm_distUnif). It must be present **if and only if** the input argument `getCDF` is missing.)
    !>  \param[in]      weight2     :   The input vector of the same size as the size of `sample2` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, <br>
    !>                                      <li>    type `real` of the same kind as that of the output `disKolm` of type `real`, <br>
    !>                                  </ol>
    !>                                  representing the weights of the corresponding elements of the input `sample2`.<br>
    !>                                  (**optional**. default = `[(1, i = 1, size(sample2))]`. It can be present only if the input argument `wieght1` and `sample2` are also present.)
    !>  \param[in]      weisum2     :   The input scalar of the same type and kind as the input `weight2`, containing the value `sum(weight2)`.<br>
    !>                                  (**optional**. It must be present if and only if the input argument `weight2` is also present.)
    !>  \param          getCDF      :   The `external` user-specified function that takes one input **scalar** argument of the same type and kind as the output `disKolm`.<br>
    !>                                  It returns a scalar of same type and kind as the output `disKolm`, representing the Cumulative Distribution Function (CDF) of the
    !>                                  continuous distribution whose Kolmogorov distance with respect to the input `sample1` must be measured.<br>
    !>                                  The following illustrates the generic interface of `getCDF()`,
    !>                                  \code{.F90}
    !>                                      function getCDF(x) result(cdf)
    !>                                          real(RKG), intent(in) :: x
    !>                                          real(RKG) :: cdf
    !>                                      end function
    !>                                  \endcode
    !>                                  where `RKG` refers to the kind type parameter of the output `disKolm`.<br>
    !>                                  (**optional**, the default is the CDF of the [Uniform distribution](@ref pm_distUnif). It must be present **if and only if** the input argument `sample2` is missing.)
    !>  \param[in]      order       :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [ascending](@ref pm_array::ascending) or
    !>                                              an object of type [ascending_type](@ref pm_array::ascending_type),
    !>                                              implying that the input samples are both sorted in ascending order.<br>
    !>                                  </ol>
    !>                                  (**optional**. If missing, the samples will be sorted in ascending order before computing the distance.)
    !>
    !>  \return
    !>  `disKolm`                   :   The output scalar of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  representing the Kolmogorov distance of `sample1` from `sample2`.<br>
    !>
    !>  \interface{getDisKolm}
    !>  \code{.F90}
    !>
    !>      use pm_distanceKolm, only: getDisKolm, ascending
    !>
    !>      ! unweighted sample against uniform distribution.
    !>
    !>      disKolm = getDisKolm(sample1(:))
    !>      disKolm = getDisKolm(sample1(:), order)
    !>
    !>      ! weighted sample against uniform distribution.
    !>
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1)
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, order)
    !>
    !>      ! unweighted sample against custom distribution CDF.
    !>
    !>      disKolm = getDisKolm(sample1(:), getCDF)
    !>      disKolm = getDisKolm(sample1(:), getCDF, order)
    !>
    !>      ! weighted sample against custom distribution CDF.
    !>
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, getCDF)
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, getCDF, order)
    !>
    !>      ! unweighted samples.
    !>
    !>      disKolm = getDisKolm(sample1(:), sample2(:))
    !>      disKolm = getDisKolm(sample1(:), sample2(:), order)
    !>
    !>      ! one weighted sample.
    !>
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, sample2(:))
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, sample2(:), order)
    !>
    !>      ! two weighted samples.
    !>
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, sample2(:), weight2(:), weisum2)
    !>      disKolm = getDisKolm(sample1(:), weight1(:), weisum1, sample2(:), weight2(:), weisum2, order)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 <= weight1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= weight2)` must hold for the corresponding input arguments.<br>
    !>  The condition `weisum1 == sum(weight1)` must hold for the corresponding input arguments.<br>
    !>  The condition `weisum2 == sum(weight2)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample1) == size(weight1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample2) == size(weight2)` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscending(sample1) .or. .not. same_type_as(order, ascending)` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscending(sample2) .or. .not. same_type_as(order, ascending)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= sample1) .and. all(sample1 <= 1) .or. present(sample2)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= sample1) .and. all(sample1 <= 1)` must hold for the corresponding input arguments when the input arguments `getCDF` and `sample2` are both missing.<br>
    !>  The condition `0 <= getCDF(x) <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The returned distance is `0` if any of the two samples are empty.<br>
    !>
    !>  \warnpure
    !>
    !>  \example{getDisKolm}
    !>  \include{lineno} example/pm_distanceKolm/getDisKolm/main.F90
    !>  \compilef{getDisKolm}
    !>  \output{getDisKolm}
    !>  \include{lineno} example/pm_distanceKolm/getDisKolm/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceKolm](@ref test_pm_distanceKolm)<br>
    !>
    !>  \naming
    !>  \code{.F90}
    !>      getDisKolmSSD_WDD_D1_RK5()
    !>                ||| ||| || |||
    !>                ||| ||| || |||
    !>                ||| ||| || |||
    !>                ||| ||| || |The Kind of the output.
    !>                ||| ||| The rank of the input sample(s).
    !>                ||| The sample weight types: WDD => default (unweighted) / default (unweighted), WID => integer-weighted, default, WRD => real-weighted, default., WII => integer-weighted, integer-weighted, WRR => real-weighted, real-weighted.
    !>                ||The order of the input samples: D => default (unordered), O => ordered.
    !>                The entities in the distance computation, S => sample, X => default (uniform) distribution, C => custom user-specified CDF.
    !>  \endcode
    !>
    !>  \final{getDisKolm}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! SS*_WDD_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSD_WDD_D1_RK5(sample1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSD_WDD_D1_RK4(sample1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSD_WDD_D1_RK3(sample1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSD_WDD_D1_RK2(sample1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSD_WDD_D1_RK1(sample1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSA_WDD_D1_RK5(sample1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSA_WDD_D1_RK4(sample1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSA_WDD_D1_RK3(sample1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSA_WDD_D1_RK2(sample1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSA_WDD_D1_RK1(sample1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WID_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSD_WID_D1_RK5(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSD_WID_D1_RK4(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSD_WID_D1_RK3(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSD_WID_D1_RK2(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSD_WID_D1_RK1(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSA_WID_D1_RK5(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSA_WID_D1_RK4(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSA_WID_D1_RK3(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSA_WID_D1_RK2(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSA_WID_D1_RK1(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WRD_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSD_WRD_D1_RK5(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSD_WRD_D1_RK4(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSD_WRD_D1_RK3(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSD_WRD_D1_RK2(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSD_WRD_D1_RK1(sample1, weight1, weisum1, sample2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSA_WRD_D1_RK5(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSA_WRD_D1_RK4(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSA_WRD_D1_RK3(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSA_WRD_D1_RK2(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSA_WRD_D1_RK1(sample1, weight1, weisum1, sample2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WII_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSD_WII_D1_RK5(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WII_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSD_WII_D1_RK4(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WII_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSD_WII_D1_RK3(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WII_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSD_WII_D1_RK2(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WII_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSD_WII_D1_RK1(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WII_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSA_WII_D1_RK5(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WII_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSA_WII_D1_RK4(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WII_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSA_WII_D1_RK3(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WII_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSA_WII_D1_RK2(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WII_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSA_WII_D1_RK1(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WII_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WRR_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSD_WRR_D1_RK5(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSD_WRR_D1_RK4(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSD_WRR_D1_RK3(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSD_WRR_D1_RK2(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSD_WRR_D1_RK1(sample1, weight1, weisum1, sample2, weight2, weisum2) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSD_WRR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSSA_WRR_D1_RK5(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSSA_WRR_D1_RK4(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSSA_WRR_D1_RK3(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSSA_WRR_D1_RK2(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSSA_WRR_D1_RK1(sample1, weight1, weisum1, sample2, weight2, weisum2, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSSA_WRR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                :: weisum1, weisum2
        type(ascending_type), intent(in)                :: order
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SX*_WDD_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSXD_WDD_D1_RK5(sample1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSXD_WDD_D1_RK4(sample1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSXD_WDD_D1_RK3(sample1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSXD_WDD_D1_RK2(sample1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSXD_WDD_D1_RK1(sample1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSXA_WDD_D1_RK5(sample1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSXA_WDD_D1_RK4(sample1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSXA_WDD_D1_RK3(sample1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSXA_WDD_D1_RK2(sample1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSXA_WDD_D1_RK1(sample1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SX*_WID_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSXD_WID_D1_RK5(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSXD_WID_D1_RK4(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSXD_WID_D1_RK3(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSXD_WID_D1_RK2(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSXD_WID_D1_RK1(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSXA_WID_D1_RK5(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSXA_WID_D1_RK4(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSXA_WID_D1_RK3(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSXA_WID_D1_RK2(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSXA_WID_D1_RK1(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SX*_WRD_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSXD_WRD_D1_RK5(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSXD_WRD_D1_RK4(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSXD_WRD_D1_RK3(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSXD_WRD_D1_RK2(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSXD_WRD_D1_RK1(sample1, weight1, weisum1) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXD_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisKolmSXA_WRD_D1_RK5(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisKolmSXA_WRD_D1_RK4(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisKolmSXA_WRD_D1_RK3(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisKolmSXA_WRD_D1_RK2(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisKolmSXA_WRD_D1_RK1(sample1, weight1, weisum1, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSXA_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SC*_WDD_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisKolmSCD_WDD_D1_RK5(sample1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK4_ENABLED
    module function getDisKolmSCD_WDD_D1_RK4(sample1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK3_ENABLED
    module function getDisKolmSCD_WDD_D1_RK3(sample1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK2_ENABLED
    module function getDisKolmSCD_WDD_D1_RK2(sample1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK1_ENABLED
    module function getDisKolmSCD_WDD_D1_RK1(sample1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisKolmSCA_WDD_D1_RK5(sample1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK4_ENABLED
    module function getDisKolmSCA_WDD_D1_RK4(sample1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK3_ENABLED
    module function getDisKolmSCA_WDD_D1_RK3(sample1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK2_ENABLED
    module function getDisKolmSCA_WDD_D1_RK2(sample1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK1_ENABLED
    module function getDisKolmSCA_WDD_D1_RK1(sample1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SC*_WID_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisKolmSCD_WID_D1_RK5(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK4_ENABLED
    module function getDisKolmSCD_WID_D1_RK4(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK3_ENABLED
    module function getDisKolmSCD_WID_D1_RK3(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK2_ENABLED
    module function getDisKolmSCD_WID_D1_RK2(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK1_ENABLED
    module function getDisKolmSCD_WID_D1_RK1(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisKolmSCA_WID_D1_RK5(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK4_ENABLED
    module function getDisKolmSCA_WID_D1_RK4(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK3_ENABLED
    module function getDisKolmSCA_WID_D1_RK3(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK2_ENABLED
    module function getDisKolmSCA_WID_D1_RK2(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK1_ENABLED
    module function getDisKolmSCA_WID_D1_RK1(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        integer(IK)         , intent(in), contiguous    :: weight1(:)
        integer(IK)         , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SC*_WRD_D1

    interface getDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisKolmSCD_WRD_D1_RK5(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK4_ENABLED
    module function getDisKolmSCD_WRD_D1_RK4(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK3_ENABLED
    module function getDisKolmSCD_WRD_D1_RK3(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK2_ENABLED
    module function getDisKolmSCD_WRD_D1_RK2(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK1_ENABLED
    module function getDisKolmSCD_WRD_D1_RK1(sample1, weight1, weisum1, getCDF) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCD_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisKolmSCA_WRD_D1_RK5(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK4_ENABLED
    module function getDisKolmSCA_WRD_D1_RK4(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK3_ENABLED
    module function getDisKolmSCA_WRD_D1_RK3(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK2_ENABLED
    module function getDisKolmSCA_WRD_D1_RK2(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

#if RK1_ENABLED
    module function getDisKolmSCA_WRD_D1_RK1(sample1, weight1, weisum1, getCDF, order) result(disKolm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisKolmSCA_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                :: order
        real(RKG)           , intent(in), contiguous    :: sample1(:)
        real(RKG)           , intent(in), contiguous    :: weight1(:)
        real(RKG)           , intent(in)                :: weisum1
        real(RKG)                                       :: disKolm
        procedure(real(RKG))                            :: getCDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Kolmogorov distance of a `sample1` of size `nsam1` from another
    !>  sample `sample2` of size `nsam2` or the CDF of the Uniform or a custom reference distribution.<br>
    !>
    !>  \details
    !>  See [pm_distanceKolm](@ref pm_distanceKolm) for the mathematical definition of the Kolmogorov distance.
    !>
    !>  \param[out]     disKolm     :   The output scalar of the same type and kind as the input `disKolm`
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  representing the Kolmogorov distance of `sample1` from `sample2`.<br>
    !>  \param[inout]   sample1     :   The input/output vector of size `nsam1` of the same type and kind as the output `disKolm`,
    !>                                  representing the first sample whose distance form the other sample must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If the `optional` input argument `order` is set to [ascending](@ref pm_array::ascending), then `sample1` has `intent(in)`.
    !>                                      <li>    If the `optional` input argument `order` is missing, then `sample1` has `intent(inout)`.<br>
    !>                                              On output, the contents of `sample1` are sorted in ascending order.<br>
    !>                                  </ol>
    !>  \param[inout]   weight1     :   The input/output vector of the same size as the size of `sample1` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, <br>
    !>                                      <li>    type `real` of the same kind as that of the output `disKolm` of type `real`, <br>
    !>                                  </ol>
    !>                                  representing the weights of the corresponding elements of the input `sample1`.<br>
    !>                                  <ol>
    !>                                      <li>    If the `optional` input argument `order` is set to [ascending](@ref pm_array::ascending), then `weight1` has `intent(in)`.
    !>                                      <li>    If the `optional` input argument `order` is missing, then `weight1` has `intent(inout)`.<br>
    !>                                              On output, the contents of `weight1` are sorted in ascending order of `sample1`.<br>
    !>                                  </ol>
    !>                                  (**optional**. default = `[(1, i = 1, size(sample1))]`. It must be present if the input argument `weight2` is also present.)
    !>  \param[in]      weisum1     :   The input scalar of the same type and kind as the input `weight1`, containing the value `sum(weight1)`.<br>
    !>                                  (**optional**. It must be present if and only if the input argument `weight1` is also present.)
    !>  \param[in]      sample2     :   The input/output vector of size `nsam2` of the same type and kind as the input `disKolm`,
    !>                                  representing the second sample whose distance form the other sample must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If the `optional` input argument `order` is set to [ascending](@ref pm_array::ascending), then `sample2` has `intent(in)`.
    !>                                      <li>    If the `optional` input argument `order` is missing, then `sample2` has `intent(inout)`.<br>
    !>                                              On output, the contents of `sample2` are sorted in ascending order.<br>
    !>                                  </ol>
    !>                                  (**optional**, the default is the [Uniform distribution](@ref pm_distUnif). It must be present **if and only if** the input argument `getCDF` is missing.)
    !>  \param[in]      weight2     :   The input vector of the same size as the size of `sample2` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, <br>
    !>                                      <li>    type `real` of the same kind as that of the output `disKolm` of type `real`, <br>
    !>                                  </ol>
    !>                                  representing the weights of the corresponding elements of the input `sample2`.<br>
    !>                                  <ol>
    !>                                      <li>    If the `optional` input argument `order` is set to [ascending](@ref pm_array::ascending), then `weight2` has `intent(in)`.
    !>                                      <li>    If the `optional` input argument `order` is missing, then `weight2` has `intent(inout)`.<br>
    !>                                              On output, the contents of `weight2` are sorted in ascending order of `sample2`.<br>
    !>                                  </ol>
    !>                                  (**optional**. default = `[(1, i = 1, size(sample2))]`. It can be present only if the input argument `wieght1` and `sample2` are also present.)
    !>  \param[in]      weisum2     :   The input scalar of the same type and kind as the input `weight2`, containing the value `sum(weight2)`.<br>
    !>                                  (**optional**. It must be present if and only if the input argument `weight2` is also present.)
    !>  \param          getCDF      :   The `external` user-specified function that takes one input **scalar** argument of the same type and kind as the output `disKolm`.<br>
    !>                                  It returns a scalar of same type and kind as the output `disKolm`, representing the Cumulative Distribution Function (CDF) of the
    !>                                  continuous distribution whose Kolmogorov distance with respect to the input `sample1` must be measured.<br>
    !>                                  The following illustrates the generic interface of `getCDF()`,
    !>                                  \code{.F90}
    !>                                      function getCDF(x) result(cdf)
    !>                                          real(RKG), intent(in) :: x
    !>                                          real(RKG) :: cdf
    !>                                      end function
    !>                                  \endcode
    !>                                  where `RKG` refers to the kind type parameter of the output `disKolm`.<br>
    !>                                  (**optional**, the default is the CDF of the [Uniform distribution](@ref pm_distUnif). It must be present **if and only if** the input argument `sample2` is missing.)
    !>  \param[in]      order       :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [ascending](@ref pm_array::ascending) or
    !>                                              an object of type [ascending_type](@ref pm_array::ascending_type),
    !>                                              implying that the input samples are both sorted in ascending order.<br>
    !>                                  </ol>
    !>                                  (**optional**. If missing, the samples will be sorted in ascending order before computing the distance.)
    !>
    !>  \interface{setDisKolm}
    !>  \code{.F90}
    !>
    !>      use pm_distanceKolm, only: setDisKolm, ascending
    !>
    !>      ! unweighted sample against uniform distribution.
    !>
    !>      call setDisKolm(disKolm, sample1(:))
    !>      call setDisKolm(disKolm, sample1(:), order)
    !>
    !>      ! weighted sample against uniform distribution.
    !>
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1)
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, order)
    !>
    !>      ! unweighted sample against custom distribution CDF.
    !>
    !>      call setDisKolm(disKolm, sample1(:), getCDF)
    !>      call setDisKolm(disKolm, sample1(:), getCDF, order)
    !>
    !>      ! weighted sample against custom distribution CDF.
    !>
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, getCDF)
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, getCDF, order)
    !>
    !>      ! unweighted samples.
    !>
    !>      call setDisKolm(disKolm, sample1(:), sample2(:))
    !>      call setDisKolm(disKolm, sample1(:), sample2(:), order)
    !>
    !>      ! one weighted sample.
    !>
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, sample2(:))
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, sample2(:), order)
    !>
    !>      ! two weighted samples.
    !>
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, sample2(:), weight2(:), weisum2)
    !>      call setDisKolm(disKolm, sample1(:), weight1(:), weisum1, sample2(:), weight2(:), weisum2, order)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 <= weight1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= weight2)` must hold for the corresponding input arguments.<br>
    !>  The condition `weisum1 == sum(weight1)` must hold for the corresponding input arguments.<br>
    !>  The condition `weisum2 == sum(weight2)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample1) == size(weight1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample2) == size(weight2)` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscending(sample1) .or. .not. same_type_as(order, ascending)` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscending(sample2) .or. .not. same_type_as(order, ascending)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= sample1) .and. all(sample1 <= 1) .or. present(sample2)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= sample1) .and. all(sample1 <= 1)` must hold for the corresponding input arguments when the input arguments `getCDF` and `sample2` are both missing.<br>
    !>  The condition `0 <= getCDF(x) <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The returned distance is `0` if any of the two samples are empty.<br>
    !>
    !>  \warnpure
    !>
    !>  \example{setDisKolm}
    !>  \include{lineno} example/pm_distanceKolm/setDisKolm/main.F90
    !>  \compilef{setDisKolm}
    !>  \output{setDisKolm}
    !>  \include{lineno} example/pm_distanceKolm/setDisKolm/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceKolm](@ref test_pm_distanceKolm)<br>
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setDisKolmSSD_WDD_D1_RK5()
    !>                ||| ||| || |||
    !>                ||| ||| || |||
    !>                ||| ||| || |||
    !>                ||| ||| || |The Kind of the output.
    !>                ||| ||| The rank of the input sample(s).
    !>                ||| The sample weight types: WDD => default (unweighted) / default (unweighted), WID => integer-weighted, default, WRD => real-weighted, default., WII => integer-weighted, integer-weighted, WRR => real-weighted, real-weighted.
    !>                ||The order of the input samples: D => default (unordered), O => ordered.
    !>                The entities in the distance computation, S => sample, X => default (uniform) distribution, C => custom user-specified CDF.
    !>  \endcode
    !>
    !>  \final{setDisKolm}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! SS*_WDD_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSD_WDD_D1_RK5(disKolm, sample1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSD_WDD_D1_RK4(disKolm, sample1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSD_WDD_D1_RK3(disKolm, sample1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSD_WDD_D1_RK2(disKolm, sample1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSD_WDD_D1_RK1(disKolm, sample1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSA_WDD_D1_RK5(disKolm, sample1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSA_WDD_D1_RK4(disKolm, sample1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSA_WDD_D1_RK3(disKolm, sample1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSA_WDD_D1_RK2(disKolm, sample1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSA_WDD_D1_RK1(disKolm, sample1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WID_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSD_WID_D1_RK5(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSD_WID_D1_RK4(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSD_WID_D1_RK3(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSD_WID_D1_RK2(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSD_WID_D1_RK1(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSA_WID_D1_RK5(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSA_WID_D1_RK4(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSA_WID_D1_RK3(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSA_WID_D1_RK2(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSA_WID_D1_RK1(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WRD_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSD_WRD_D1_RK5(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSD_WRD_D1_RK4(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSD_WRD_D1_RK3(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSD_WRD_D1_RK2(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSD_WRD_D1_RK1(disKolm, sample1, weight1, weisum1, sample2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSA_WRD_D1_RK5(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSA_WRD_D1_RK4(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSA_WRD_D1_RK3(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSA_WRD_D1_RK2(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSA_WRD_D1_RK1(disKolm, sample1, weight1, weisum1, sample2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WII_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSD_WII_D1_RK5(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WII_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSD_WII_D1_RK4(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WII_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSD_WII_D1_RK3(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WII_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSD_WII_D1_RK2(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WII_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSD_WII_D1_RK1(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WII_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSA_WII_D1_RK5(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WII_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSA_WII_D1_RK4(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WII_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSA_WII_D1_RK3(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WII_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSA_WII_D1_RK2(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WII_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSA_WII_D1_RK1(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WII_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:), weight2(:)
        integer(IK)         , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SS*_WRR_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSD_WRR_D1_RK5(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSD_WRR_D1_RK4(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSD_WRR_D1_RK3(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSD_WRR_D1_RK2(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSD_WRR_D1_RK1(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSD_WRR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSSA_WRR_D1_RK5(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSSA_WRR_D1_RK4(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSSA_WRR_D1_RK3(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSSA_WRR_D1_RK2(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSSA_WRR_D1_RK1(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSSA_WRR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:), sample2(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:), weight2(:)
        real(RKG)           , intent(in)                    :: weisum1, weisum2
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SX*_WDD_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSXD_WDD_D1_RK5(disKolm, sample1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSXD_WDD_D1_RK4(disKolm, sample1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSXD_WDD_D1_RK3(disKolm, sample1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSXD_WDD_D1_RK2(disKolm, sample1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSXD_WDD_D1_RK1(disKolm, sample1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSXA_WDD_D1_RK5(disKolm, sample1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSXA_WDD_D1_RK4(disKolm, sample1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSXA_WDD_D1_RK3(disKolm, sample1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSXA_WDD_D1_RK2(disKolm, sample1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSXA_WDD_D1_RK1(disKolm, sample1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(ascending_type), intent(in)                    :: order
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SX*_WID_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSXD_WID_D1_RK5(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSXD_WID_D1_RK4(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSXD_WID_D1_RK3(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSXD_WID_D1_RK2(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSXD_WID_D1_RK1(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSXA_WID_D1_RK5(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSXA_WID_D1_RK4(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSXA_WID_D1_RK3(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSXA_WID_D1_RK2(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSXA_WID_D1_RK1(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SX*_WRD_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSXD_WRD_D1_RK5(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSXD_WRD_D1_RK4(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSXD_WRD_D1_RK3(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSXD_WRD_D1_RK2(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSXD_WRD_D1_RK1(disKolm, sample1, weight1, weisum1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXD_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisKolmSXA_WRD_D1_RK5(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisKolmSXA_WRD_D1_RK4(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisKolmSXA_WRD_D1_RK3(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisKolmSXA_WRD_D1_RK2(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisKolmSXA_WRD_D1_RK1(disKolm, sample1, weight1, weisum1, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSXA_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SC*_WDD_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setDisKolmSCD_WDD_D1_RK5(disKolm, sample1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setDisKolmSCD_WDD_D1_RK4(disKolm, sample1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setDisKolmSCD_WDD_D1_RK3(disKolm, sample1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setDisKolmSCD_WDD_D1_RK2(disKolm, sample1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setDisKolmSCD_WDD_D1_RK1(disKolm, sample1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setDisKolmSCA_WDD_D1_RK5(disKolm, sample1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setDisKolmSCA_WDD_D1_RK4(disKolm, sample1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setDisKolmSCA_WDD_D1_RK3(disKolm, sample1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setDisKolmSCA_WDD_D1_RK2(disKolm, sample1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setDisKolmSCA_WDD_D1_RK1(disKolm, sample1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SC*_WID_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setDisKolmSCD_WID_D1_RK5(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setDisKolmSCD_WID_D1_RK4(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setDisKolmSCD_WID_D1_RK3(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setDisKolmSCD_WID_D1_RK2(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setDisKolmSCD_WID_D1_RK1(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        integer(IK)         , intent(inout) , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setDisKolmSCA_WID_D1_RK5(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setDisKolmSCA_WID_D1_RK4(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setDisKolmSCA_WID_D1_RK3(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setDisKolmSCA_WID_D1_RK2(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setDisKolmSCA_WID_D1_RK1(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        integer(IK)         , intent(in)    , contiguous    :: weight1(:)
        integer(IK)         , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SC*_WRD_D1

    interface setDisKolm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setDisKolmSCD_WRD_D1_RK5(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setDisKolmSCD_WRD_D1_RK4(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setDisKolmSCD_WRD_D1_RK3(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setDisKolmSCD_WRD_D1_RK2(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setDisKolmSCD_WRD_D1_RK1(disKolm, sample1, weight1, weisum1, getCDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCD_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample1(:)
        real(RKG)           , intent(inout) , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setDisKolmSCA_WRD_D1_RK5(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WRD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setDisKolmSCA_WRD_D1_RK4(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WRD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setDisKolmSCA_WRD_D1_RK3(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WRD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setDisKolmSCA_WRD_D1_RK2(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WRD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setDisKolmSCA_WRD_D1_RK1(disKolm, sample1, weight1, weisum1, getCDF, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisKolmSCA_WRD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: sample1(:)
        real(RKG)           , intent(in)    , contiguous    :: weight1(:)
        real(RKG)           , intent(in)                    :: weisum1
        real(RKG)           , intent(out)                   :: disKolm
        type(ascending_type), intent(in)                    :: order
        procedure(real(RKG))                                :: getCDF
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distanceKolm ! LCOV_EXCL_LINE
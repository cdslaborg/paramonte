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
!>  This module contains procedures and data types for computing sample quantile.<br>
!>
!>  \details
!>  Quantiles are cut points dividing the range of a probability distribution into continuous
!>  intervals with equal probabilities, or dividing the observations in a sample in the same way.<br>
!>  There is one fewer quantile than the number of groups created.<br>
!>  Common quantiles have special names, such as quartiles (four groups), deciles (ten groups), and percentiles (100 groups).<br>
!>  The groups created are termed halves, thirds, quarters, etc., though sometimes the terms for the quantile are used for the groups created, rather than for the cut points.<br>
!>
!>  Q-quantiles are values that partition a finite set of values into \f$q\f$ subsets of (nearly) equal sizes.<br>
!>  There are \f$q − 1\f$ partitions of the \f$q\f$-quantiles, one for each integer \f$k\f$ satisfying \f$0 < k < q\f$.<br>
!>  In some cases the value of a quantile may not be uniquely determined, as can be the case for the median (\f$2\f$-quantile) of a uniform probability distribution on a set of even size.<br>
!>  Quantiles can also be applied to continuous distributions, providing a way to generalize rank statistics to continuous variables.<br>
!>  When the cumulative distribution function of a random variable is known, the \f$q\f$-quantiles are the application of the quantile function
!>  (the inverse function of the cumulative distribution function) to the values \f$\{1/q, 2/q, …, (q − 1)/q\}\f$.<br>
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
!>  [pm_polation](@ref pm_polation)<br>
!>
!>  \test
!>  [test_pm_sampleQuan](@ref test_pm_sampleQuan)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleQuan

    use pm_polation, only: setInterp
    use pm_polation, only: neimean, neimean_type
    use pm_polation, only: neinear, neinear_type
    use pm_polation, only: neinext, neinext_type
    use pm_polation, only: neiprev, neiprev_type
    use pm_polation, only: piwilin, piwilin_type
    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleQuan"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the approximate sample quantile for the given `method` at the specified probabilities.
    !>
    !>  \param[in]  method  :   The input scalar constant that can be,
    !>                          <ol>
    !>                              <li>    The scalar constant [neimean](@ref pm_polation::neimean) or a scalar object of type [neimean_type](@ref pm_polation::neimean_type)
    !>                                      implying the use of the average of the `sample` values of the two nearest neighbors of the input `prob` smaller and larger than it as the output `quan`.<br>
    !>                              <li>    The scalar constant [neinear](@ref pm_polation::neinear) or a scalar object of type [neinear_type](@ref pm_polation::neinear_type)
    !>                                      implying the use of the average of the `sample` value of the nearest neighbor of the input `prob` as the output `quan`.<br>
    !>                                      Note that the nearest neighbor in this case is measured by actual Euclidean distances of neighbors to the input `prob`.<br>
    !>                              <li>    The scalar constant [neiprev](@ref pm_polation::neiprev) or a scalar object of type [neiprev_type](@ref pm_polation::neiprev_type)
    !>                                      implying the use of the `sample` value of the largest abscissa in the input `prob` smaller than the input `prob` as the output `quan`.<br>
    !>                              <li>    The scalar constant [neinext](@ref pm_polation::neinext) or a scalar object of type [neinext_type](@ref pm_polation::neinext_type)
    !>                                      implying the use of the `sample` value of the smallest abscissa in the input `prob` larger than the input `prob` as the output `quan`.<br>
    !>                              <li>    The scalar constant [piwilin](@ref pm_polation::piwilin) or a scalar object of type [piwilin_type](@ref pm_polation::piwilin_type)
    !>                                      implying the use of the **linear interpolation** of the `sample` values of the two `prob` points that bracket `prob` as the output `quan`.<br>
    !>                                      The linear interpolation implemented in this constructor is based on the Lagrange classical formula for linear interpolation.<br>
    !>                                      Suppose an input query point \f$x\f$ falls between two nodes \f$x_i\f$ and \f$x_{i+1}\f$ with the corresponding function values
    !>                                      \f$y_i\f$ and \f$y_{i+1}\f$ and we wish to estimate the corresponding interpolated value \f$y(x)\f$, which can be computed as,
    !>                                      \f{equation*}{
    !>                                           y(x) = \frac {x - x_{i+1}} {x_i - x_{i+1}} y_i + \frac {x - x_{i}} {x_{i+1} - x_{i}} y_{i+1} ~.
    !>                                      \f}
    !>                          </ol>
    !>                          Either [piwilin](@ref pm_polation::piwilin) or [neinear](@ref pm_polation::neinear) for larger sample sizes can be reasonable choices.<br>
    !>  \param[in]  prob    :   The input scalar or `contiguous` vector of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the (set of) probability(s) for which the corresponding quantiles must be computed.<br>
    !>  \param[in]  sample  :   The input `contiguous` vector of shape `(1:nsam)` or matrix of shape `(1:ndim, 1:nsam)` of the same type, kind, and size as `prob`,
    !>                          containing the sample based upon which the quantiles must be computed.
    !>  \param[in]  dim     :   The input scalar of type `integer` of default kind \IK,
    !>                          whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                          <ol>
    !>                              <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(1:nsam, 1:ndim)`.<br>
    !>                              <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(1:ndim, 1:nsam)`.<br>
    !>                          </ol>
    !>                          (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]  weight  :   The input `contiguous` vector of length `nsam` of,
    !>                          <ol>
    !>                              <li>    type `integer` of default kind \IK, or
    !>                              <li>    type `real` of the same kind as the kind of `prob`,
    !>                          </ol>
    !>                          containing the corresponding weights of individual `nsam` observations in the target sample.<br>
    !>                          (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled).)
    !>  \param[in]  weisum  :   The input scalar of the same type and kind as `weight` containing the quantity `sum(weight)`.<br>
    !>                          (**optional**. It **can** be present **if and only if** the input argument `weight` is also present.)
    !>
    !>  \return
    !>
    !>  `quan`              :   The output object of the same type and kind as `prob`,
    !>                          containing the sample quantiles corresponding to the input `prob`.
    !>                          <ol>
    !>                              <li>    If `sample` is a vector and `prob` is a scalar, then `quan` is a scalar.
    !>                              <li>    If `sample` is a vector and `prob` is a vector, then `quan` is a vector of size `size(prob)`.
    !>                              <li>    If `sample` is a matrix and `prob` is a scalar, then `quan` is a vector of size `size(sample, 3 - dim)`.
    !>                              <li>    If `sample` is a matrix and `prob` is a vector, then `quan` is a matrix of shape `[size(prob), size(sample, 3 - dim)]`.
    !>                          </ol>
    !>
    !>  \interface{getQuan}
    !>  \code{.F90}
    !>
    !>      use pm_sampleQuan, only: getQuan
    !>
    !>      ! 1D sample
    !>
    !>      quan = getQuan(method, prob, sample(1:nsam))
    !>      quan = getQuan(method, prob, sample(1:nsam), weight(1:nsam), weisum = weisum)
    !>
    !>      quan(1:nprob) = getQuan(method, prob(1:nprob), sample(1:nsam))
    !>      quan(1:nprob) = getQuan(method, prob(1:nprob), sample(1:nsam), weight(1:nsam), weisum = weisum)
    !>
    !>      ! 2D sample
    !>
    !>      quan(1:size(sample, 3 - dim)) = getQuan(method, prob, sample(:,:), dim)
    !>      quan(1:size(sample, 3 - dim)) = getQuan(method, prob, sample(:,:), dim, weight(1:nsam), weisum = weisum)
    !>
    !>      quan(1:nprob, 1:size(sample, 3 - dim)) = getQuan(method, prob(1:nprob), sample(:,:), dim)
    !>      quan(1:nprob, 1:size(sample, 3 - dim)) = getQuan(method, prob(1:nprob), sample(:,:), dim, weight(1:nsam), weisum = weisum)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(weight) == size(sample, dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `all([0 <= prob .and. prob <= 1])` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < dim .and. dim < 3` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= weight)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \devnote
    !>  A subroutine equivalent of this functional interface was deemed unnecessary.<br>
    !>  Internally, the functional interface requires two internal runtime `allocatable` arrays to sort the input sample and store the sample ECDF.<br>
    !>  For intensive repeated calculations, these `allocation` actions and copies can be avoided by explicitly implementing this generic interface as the following for a 1D sample.<br>
    !>  \code{.F90}
    !>      use pm_arraySort, only: setSorted
    !>      use pm_sampleECDF, only: setECDF
    !>      use pm_polation, only: setExtrap
    !>      type_of(sample) :: ecdf(size(sample))
    !>      call setECDF(ecdf, weight, weisum)
    !>      call setSorted(sample)
    !>      call setExtrap(method, ecdf, sample, prob, quan)
    !>  \endcode
    !>
    !>  \see
    !>  [pm_polation](@ref pm_polation)<br>
    !>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
    !>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
    !>  [pm_arraySort](@ref pm_arraySort)<br>
    !>
    !>  \example{getQuan}
    !>  \include{lineno} example/pm_sampleQuan/getQuan/main.F90
    !>  \compilef{getQuan}
    !>  \output{getQuan}
    !>  \include{lineno} example/pm_sampleQuan/getQuan/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleQuan](@ref test_pm_sampleQuan)
    !>
    !>  \final{getQuan}
    !>
    !>  \author
    !>  Amir Shahmoradi, Thursday 2:45 AM, August 19, 2021, Dallas, TX<br>

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ND1 WNO piwilin

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WNO neimean

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WNO neinear

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WNO neinext

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WNO neiprev

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WNO_RK5(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WNO_RK4(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WNO_RK3(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WNO_RK2(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WNO_RK1(method, prob, sample) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ND1 WTI piwilin

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTI neimean

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTI neinear

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTI neinext

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTI neiprev

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTI_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTI_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTI_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTI_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTI_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ND1 WTR piwilin

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND1_QD0_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND1_QD1_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND1_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTR neimean

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND1_QD0_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND1_QD1_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND1_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTR neinear

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND1_QD0_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND1_QD1_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND1_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTR neinext

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND1_QD0_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND1_QD1_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND1_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND1 WTR neiprev

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND1_QD0_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)                                           :: quan
        real(TKG)           , intent(in)                    :: prob
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTR_RK5(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTR_RK4(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTR_RK3(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTR_RK2(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND1_QD1_WTR_RK1(method, prob, sample, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND1_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: sample(:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ND2 WNO piwilin

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WNO neimean

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WNO neinear

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WNO neinext

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WNO neiprev

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WNO_RK5(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WNO_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WNO_RK4(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WNO_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WNO_RK3(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WNO_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WNO_RK2(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WNO_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WNO_RK1(method, prob, sample, dim) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WNO_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ND2 WTI piwilin

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTI neimean

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTI neinear

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTI neinext

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTI neiprev

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTI_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTI_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTI_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTI_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTI_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTI_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTI_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTI_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTI_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTI_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        integer(IK)         , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ND2 WTR piwilin

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND2_QD0_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPWLN_ND2_QD1_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPWLN_ND2_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTR neimean

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND2_QD0_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanMEAN_ND2_QD1_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanMEAN_ND2_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTR neinear

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND2_QD0_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEAR_ND2_QD1_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEAR_ND2_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTR neinext

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND2_QD0_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanNEXT_ND2_QD1_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanNEXT_ND2_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ND2 WTR neiprev

    interface getQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND2_QD0_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD0_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        real(TKG)           , intent(in)                    :: prob
        real(TKG)                                           :: quan(size(sample, 3 - dim, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTR_RK5(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTR_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTR_RK4(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTR_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTR_RK3(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTR_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTR_RK2(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTR_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getQuanPREV_ND2_QD1_WTR_RK1(method, prob, sample, dim, weight, weisum) result(quan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuanPREV_ND2_QD1_WTR_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(TKG)           , intent(in)    , contiguous    :: sample(:,:)
        real(TKG)           , intent(in)    , contiguous    :: prob(:)
        real(TKG)                                           :: quan(size(prob, 1, IK), size(sample, 3 - dim, IK))
        real(TKG)           , intent(in)    , contiguous    :: weight(:)
        real(TKG)           , intent(in)    , optional      :: weisum
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleQuan
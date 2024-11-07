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
!>  This module contains classes and procedures for computing the first moment (i.e., the statistical mean) of random weighted samples.
!>
!>  \details
!>  The mean of a weighted sample of \f$N\f$ data points is computed by the following equation,
!>  \f{equation}{
!>      \mu = \frac{\sum_{i = 1}^{N} w_i x_i}{\sum_{i = 1}^{N} w_i}
!>  \f}
!>  where \f$w_i\f$ represents the weight of the \f$i\f$th sample.
!>
!>  Mean updating
!>  -------------
!>
!>  Suppose the mean of an initial (potentially weighted) sample \f$x_A\f$ of size \f$N_A\f$ is computed to be \f$\mu_A\f$.<br>
!>  Another (potentially weighted) sample \f$x_B\f$ of size \f$N_B\f$ is subsequently obtained with a different number observations and mean \f$\mu_B\f$.<br>
!>  The mean of the two samples combined can be expressed in terms of the originally computed means,
!>  \f{equation}{
!>      \large
!>      \mu = \frac
!>      {
!>          w_A \sum_{i = 1}^{N_A} w_{\up{A,i}} x_{\up{A,i}} +
!>          w_B \sum_{i = 1}^{N_B} w_{\up{B,i}} x_{\up{B,i}}
!>      }{
!>          w_A + w_B
!>      }
!>  \f}
!>  where \f$\large w_A = \sum_{i = 1}^{N_A} w_{\up{A,i}}\f$ and \f$\large w_B = \sum_{i = 1}^{N_B} w_{\up{B,i}}\f$
!>  are sums of the weights of the corresponding samples.<br>
!>  For equally-weighted samples, the corresponding weights \f$w_{\up{A,i}}\f$ or \f$w_{\up{B,i}}\f$ or both are all unity
!>  such that \f$N_A = w_A\f$ or \f$N_B = w_B\f$ or both holds.<br>
!>
!>  \devnote
!>  While it is tempting to extend the generic interfaces of this module to `weight` arguments of type `integer` or `real` of various kinds,
!>  such extensions do not add any benefits beyond making the interface more flexible for the end user.<br>
!>  But such extensions would certainly make the maintenance and future extensions of this interface difficult and complex.<br>
!>  According to the coercion rules of the Fortran standard, if an `integer` is multiplied with a `real`,
!>  the integer value must be first converted to `real` of the same kind as the real value, then multiplied.<br>
!>  Furthermore, the floating-point multiplication tends to be faster than integer multiplication on most modern architecture.<br>
!>  The following list compares the cost and latencies of some of basic operations involving integers and real numbers.<br>
!>  <ol>
!>      <li>    Central Processing Unit (CPU):
!>              <ol>
!>                  <li>    Integer add: 1 cycle
!>                  <li>    32-bit integer multiply: 10 cycles
!>                  <li>    64-bit integer multiply: 20 cycles
!>                  <li>    32-bit integer divide: 69 cycles
!>                  <li>    64-bit integer divide: 133 cycles
!>              </ol>
!>      <li>    On-chip Floating Point Unit (FPU):
!>              <ol>
!>                  <li>    Floating point add: 4 cycles
!>                  <li>    Floating point multiply: 7 cycles
!>                  <li>    Double precision multiply: 8 cycles
!>                  <li>    Floating point divide: 23 cycles
!>                  <li>    Double precision divide: 36 cycles
!>              </ol>
!>  </ol>
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
!>  [Intel Fortran Forum - Integer VS fp performance](https://community.intel.com/t5/Intel-Fortran-Compiler/Integer-VS-fp-performance/m-p/756639)<br>
!>  [Colorado State University tips on Fortran performance](http://hogback.atmos.colostate.edu/rr/old/tidbits/intel/macintel/doc_files/source/extfile/optaps_for/fortran/optaps_prg_runt_f.htm)<br>
!>
!>  \test
!>  [test_pm_sampleMean](@ref test_pm_sampleMean)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleMean

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleMean"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (weighted) mean of an input `sample` of `nsam`
    !>  observations with `ndim = 1 or 2` attributes, optionally weighted by the input `weight`.<br>
    !>
    !>  \param[in]  sample  :   The `contiguous` array of shape `(nsam)`, `(ndim, nsam)` or `(nsam, ndim)` of,
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the sample whose `mean` is to be computed.<br>
    !>  \param[in]  dim     :   The input scalar `integer` of default kind \IK representing the dimension (`1` or `2`) of the input `sample` along which the mean must be computed.<br>
    !>                          <ol>
    !>                              <li>    If `dim = 1`, the input `sample` of rank `2` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                              <li>    If `dim = 2`, the input `sample` of rank `2` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                          </ol>
    !>                          The input `dim` must always be `1` or missing for an input `sample` of rank `1`.<br>
    !>                          (**optional**. If missing, the mean of the whole input `sample` is computed.)
    !>  \param[in]  weight  :   The `contiguous` vector of length `nsam` of
    !>                          <ol>
    !>                              <li>    type `integer` of default kind \IK, or
    !>                              <li>    type `real` of the same kind as that of the output `mean`,
    !>                          </ol>
    !>                          containing the corresponding weights of the data points in the input `sample`.<br>
    !>                          (**optional**, default = a vector of ones.)
    !>
    !>  \return
    !>  `mean`              :   The output object of the same type and kind as the input `sample`.<br>
    !>                          <ol>
    !>                              <li>    When the input `sample` has the shape `(nsam)`, or `(:,:)` and `dim` is missing, the output mean is a scalar.<br>
    !>                              <li>    When the input `sample` has the shape `(nsam, ndim)` or `(ndim, nsam)`, the output mean is a vector of length `ndim`.<br>
    !>                          </ol>
    !>
    !>  \interface{getMean}
    !>  \code{.F90}
    !>
    !>      use pm_sampleMean, only: getMean
    !>
    !>      ! 1D sample.
    !>
    !>      mean = getMean(sample(1 : nsam))
    !>      mean = getMean(sample(1 : nsam), weight(1 : nsam))
    !>
    !>      mean = getMean(sample(1 : nsam), dim)
    !>      mean = getMean(sample(1 : nsam), dim, weight(1 : nsam))
    !>
    !>      ! 2D sample.
    !>
    !>      mean(1 : ndim) = getMean(sample(:,:))
    !>      mean(1 : ndim) = getMean(sample(:,:), weight(1 : nsam))
    !>
    !>      mean(1 : ndim) = getMean(sample(:,:), dim)
    !>      mean(1 : ndim) = getMean(sample(:,:), dim, weight(1 : nsam))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, dim) == size(weight, 1)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  If the input sample is to be an array of type `integer`, simply convert the sample to
    !>  an array of type `real` of the desired kind for the output `real` mean of the sample.<br>
    !>  There is no point in accepting an input sample of type `integer` since it will be inevitably
    !>  converted to an array of type `real` within the procedure to avoid potential integer overflow.<br>
    !>  Furthermore, an `sample` of type `integer` creates ambiguity about the `kind` of the `real`-valued returned mean by the procedure.<br>
    !>  It is therefore sensible to offer interfaces for only weights of type `integer` of default kind and `real` of the same kind as the sample.
    !>
    !>  \note
    !>  Note that the mean of any one or two-dimensional sample can be simply computed via the Fortran intrinsic routine `sum()`:
    !>  \code{.F90}
    !>      integer                         :: i
    !>      integer , parameter             :: NDIM = 3_IK
    !>      integer , parameter             :: NSAM = 1000_IK
    !>      real    , parameter             :: sample(NDIM,NSAM) = reshape([( real(i,RK), i = 1, NSAM )], shape = shape(sample))
    !>      real    , allocatable   :: mean(:)
    !>      mean = sum(sample, dim = 1) / size(transpose(sample), dim = 1)  ! assuming the first dimension represents the observations.
    !>      mean = sum(sample, dim = 2) / size(sample, dim = 2)             ! assuming the second dimension represents the observations.
    !>  \endcode
    !>
    !>  \note
    !>  The mean of a whole multidimensional array can be obtained by either
    !>  <ol>
    !>      <li>    reshaping the array to a vector form and passing it to this procedure, or
    !>      <li>    mapping the array to a 1-dimensional pointer array of the same size as the `ndim` dimensional array.
    !>  </ol>
    !>  See the examples below.<br>
    !>
    !>  \devnote
    !>  An `XY` input sample interface is impossible due to ambiguity with existing interfaces.<br>
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [setVarMean](@ref pm_sampleVar::setVarMean)<br>
    !>
    !>  \example{getMean}
    !>  \include{lineno} example/pm_sampleMean/getMean/main.F90
    !>  \compilef{getMean}
    !>  \output{getMean}
    !>  \include{lineno} example/pm_sampleMean/getMean/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleMean](@ref test_pm_sampleMean)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! ALL D1

    interface getMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanALL_WNO_D1_CK5(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanALL_WNO_D1_CK4(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanALL_WNO_D1_CK3(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanALL_WNO_D1_CK2(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanALL_WNO_D1_CK1(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanALL_WNO_D1_RK5(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanALL_WNO_D1_RK4(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanALL_WNO_D1_RK3(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanALL_WNO_D1_RK2(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanALL_WNO_D1_RK1(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanALL_WTI_D1_CK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanALL_WTI_D1_CK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanALL_WTI_D1_CK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanALL_WTI_D1_CK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanALL_WTI_D1_CK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanALL_WTI_D1_RK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanALL_WTI_D1_RK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanALL_WTI_D1_RK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanALL_WTI_D1_RK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanALL_WTI_D1_RK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanALL_WTR_D1_CK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanALL_WTR_D1_CK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanALL_WTR_D1_CK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanALL_WTR_D1_CK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanALL_WTR_D1_CK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanALL_WTR_D1_RK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanALL_WTR_D1_RK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanALL_WTR_D1_RK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanALL_WTR_D1_RK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanALL_WTR_D1_RK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ALL D2

    interface getMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanALL_WNO_D2_CK5(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanALL_WNO_D2_CK4(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanALL_WNO_D2_CK3(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanALL_WNO_D2_CK2(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanALL_WNO_D2_CK1(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanALL_WNO_D2_RK5(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanALL_WNO_D2_RK4(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanALL_WNO_D2_RK3(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanALL_WNO_D2_RK2(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanALL_WNO_D2_RK1(sample) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanALL_WTI_D2_CK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanALL_WTI_D2_CK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanALL_WTI_D2_CK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanALL_WTI_D2_CK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanALL_WTI_D2_CK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanALL_WTI_D2_RK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanALL_WTI_D2_RK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanALL_WTI_D2_RK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanALL_WTI_D2_RK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanALL_WTI_D2_RK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanALL_WTR_D2_CK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanALL_WTR_D2_CK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanALL_WTR_D2_CK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanALL_WTR_D2_CK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanALL_WTR_D2_CK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanALL_WTR_D2_RK5(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanALL_WTR_D2_RK4(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanALL_WTR_D2_RK3(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanALL_WTR_D2_RK2(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanALL_WTR_D2_RK1(sample, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanALL_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DIM D1

    interface getMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanDIM_WNO_D1_CK5(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanDIM_WNO_D1_CK4(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanDIM_WNO_D1_CK3(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanDIM_WNO_D1_CK2(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanDIM_WNO_D1_CK1(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanDIM_WNO_D1_RK5(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanDIM_WNO_D1_RK4(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanDIM_WNO_D1_RK3(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanDIM_WNO_D1_RK2(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanDIM_WNO_D1_RK1(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanDIM_WTI_D1_CK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanDIM_WTI_D1_CK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanDIM_WTI_D1_CK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanDIM_WTI_D1_CK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanDIM_WTI_D1_CK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanDIM_WTI_D1_RK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanDIM_WTI_D1_RK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanDIM_WTI_D1_RK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanDIM_WTI_D1_RK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanDIM_WTI_D1_RK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanDIM_WTR_D1_CK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanDIM_WTR_D1_CK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanDIM_WTR_D1_CK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanDIM_WTR_D1_CK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanDIM_WTR_D1_CK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanDIM_WTR_D1_RK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanDIM_WTR_D1_RK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanDIM_WTR_D1_RK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanDIM_WTR_D1_RK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanDIM_WTR_D1_RK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DIM D2

    interface getMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanDIM_WNO_D2_CK5(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanDIM_WNO_D2_CK4(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanDIM_WNO_D2_CK3(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanDIM_WNO_D2_CK2(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanDIM_WNO_D2_CK1(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanDIM_WNO_D2_RK5(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanDIM_WNO_D2_RK4(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanDIM_WNO_D2_RK3(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanDIM_WNO_D2_RK2(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanDIM_WNO_D2_RK1(sample, dim) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanDIM_WTI_D2_CK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanDIM_WTI_D2_CK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanDIM_WTI_D2_CK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanDIM_WTI_D2_CK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanDIM_WTI_D2_CK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanDIM_WTI_D2_RK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanDIM_WTI_D2_RK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanDIM_WTI_D2_RK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanDIM_WTI_D2_RK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanDIM_WTI_D2_RK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanDIM_WTR_D2_CK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanDIM_WTR_D2_CK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanDIM_WTR_D2_CK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanDIM_WTR_D2_CK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanDIM_WTR_D2_CK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)                                                :: mean(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanDIM_WTR_D2_RK5(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanDIM_WTR_D2_RK4(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanDIM_WTR_D2_RK3(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanDIM_WTR_D2_RK2(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanDIM_WTR_D2_RK1(sample, dim, weight) result(mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanDIM_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)                                                   :: mean(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the (weighted) mean of a pair of time series or of an input `sample` of `nsam`
    !>  observations with `ndim = 1 or 2` attributes, optionally weighted by the input `weight`,
    !>  optionally also `sum(weight)` and optionally, `sum(weight**2)`.<br>
    !>
    !>  \details
    !>  This generic interface is developed to specifically improve the performance of other procedures
    !>  where `sum(weight)` or `sum(weight**2)` are also needed in addition to the mean of the weighted sample.<br>
    !>  For example, this situation occurs frequently in the computation of the variance or covariance matrix of a weighted sample.<br>
    !>
    !>  \param[out] mean    :   The output object of,
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the data mean.<br>
    !>                          <ol>
    !>                              <li>    When the input `sample` has the shape `(nsam)`, the output mean must be a scalar.<br>
    !>                              <li>    When the input `x` and `y` are present, the output `mean` must be of shape `mean(1:2)`.<br>
    !>                              <li>    When the input `sample` has the shape `(nsam, ndim)` or `(ndim, nsam)`, the output mean must be a vector of length `ndim`.<br>
    !>                          </ol>
    !>  \param[in]  x       :   The `contiguous` array of shape `(nsam)` of the same type and kind as the output `mean`,
    !>                          containing the first of the two data series whose `mean` is to be computed.<br>
    !>                          (**optional**. It must be present *if and only if** the input argument `sample` is missing and `y` is present.)
    !>  \param[in]  y       :   The `contiguous` of the same type and kind and shape as the input `x`,
    !>                          containing the second of the two data series whose `mean` is to be computed.<br>
    !>                          (**optional**. It must be present *if and only if** the input argument `sample` is missing and `x` is present.)
    !>  \param[in]  sample  :   The `contiguous` array of shape `(nsam)`, `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `mean`,
    !>                          containing the sample whose `mean` is to be computed.<br>
    !>                          (**optional**. It must be present *if and only if** the input arguments `x` and `y` are missing.)
    !>  \param[in]  dim     :   The input scalar `integer` of default kind \IK representing the dimension (`1` or `2`) of the input `sample` along which the mean must be computed.<br>
    !>                          <ol>
    !>                              <li>    If `dim = 1`, the input `sample` of rank `2` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                              <li>    If `dim = 2`, the input `sample` of rank `2` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                          </ol>
    !>                          The input `dim` must always be `1` or missing for an input `sample` of rank `1`.<br>
    !>                          (**optional**. It can be present only if `sample` is also present. If missing, the mean of the whole input `sample` is computed.)
    !>  \param[in]  weight  :   The `contiguous` vector of length `nsam` of
    !>                          <ol>
    !>                              <li>    type `integer` of default kind \IK, or
    !>                              <li>    type `real` of the same kind as that of the output `mean`,
    !>                          </ol>
    !>                          containing the corresponding weights of the data points in the input `sample` or the input pair `x` and `y`.<br>
    !>  \param[out] weisum  :   The output scalar of the same type and kind as the input `weight`.<br>
    !>                          On output, it will contain `sum(weight)`.<br>
    !>                          This quantity is frequently needed in computing the weighted sample variance.<br>
    !>                          (**optional**. It must be present **if and only if** the input argument `weight` is also present.)
    !>
    !>  \interface{setMean}
    !>  \code{.F90}
    !>
    !>      use pm_sampleMean, only: setMean
    !>
    !>      ! XY sample.
    !>
    !>      call setMean(mean(1:2), x(1 : nsam), y(1 : nsam))
    !>      call setMean(mean(1:2), x(1 : nsam), y(1 : nsam), weight(1 : nsam), weisum)
    !>
    !>      ! 1D sample.
    !>
    !>      call setMean(mean, sample(1 : nsam))
    !>      call setMean(mean, sample(1 : nsam), weight(1 : nsam), weisum)
    !>
    !>      call setMean(mean, sample(1 : nsam), dim)
    !>      call setMean(mean, sample(1 : nsam), dim, weight(1 : nsam), weisum)
    !>
    !>      ! 2D sample.
    !>
    !>      call setMean(mean(1 : ndim), sample(:,:))
    !>      call setMean(mean(1 : ndim), sample(:,:), weight(1 : nsam), weisum)
    !>
    !>      call setMean(mean(1 : ndim), sample(:,:), dim)
    !>      call setMean(mean(1 : ndim), sample(:,:), dim, weight(1 : nsam), weisum)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, dim) == size(weight, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `(dim == 1 .and. size(mean, 1) == size(sample, 2)) .or. (dim == 2 .and. size(mean, 1) == size(sample, 1))` must hold for the corresponding input arguments.<br>
    !>  The condition `size(x) == size(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(x) == size(y)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  If the input sample is to be an array of type `integer`, simply convert the sample to
    !>  an array of type `real` of the desired kind for the output `real` mean of the sample.<br>
    !>  There is no point in accepting an input sample of type `integer` since it will be inevitably
    !>  converted to an array of type `real` within the procedure to avoid potential integer overflow.<br>
    !>  Furthermore, an `sample` of type `integer` creates ambiguity about the `kind` of the `real`-valued returned mean by the procedure.<br>
    !>  See the notes in the description of the [pm_sampleMean](@ref pm_sampleMean).<br>
    !>
    !>  \note
    !>  Note that the mean of any one or two-dimensional sample can be simply computed via the Fortran intrinsic routine `sum()`:
    !>  \code{.F90}
    !>      integer                         :: i
    !>      integer , parameter             :: NDIM = 3_IK
    !>      integer , parameter             :: NSAM = 1000_IK
    !>      real    , parameter             :: sample(NDIM,NSAM) = reshape([( real(i,RK), i = 1, NSAM )], shape = shape(sample))
    !>      real    , allocatable   :: mean(:)
    !>      mean = sum(sample, dim = 1) / size(transpose(sample), dim = 1)  ! assuming the first dimension represents the observations.
    !>      mean = sum(sample, dim = 2) / size(sample, dim = 2)             ! assuming the second dimension represents the observations.
    !>  \endcode
    !>
    !>  \note
    !>  The mean of a whole multidimensional array can be obtained by either,
    !>  <ol>
    !>      <li>    reshaping the array to a vector form and passing it to this procedure, or
    !>      <li>    mapping the array to a 1-dimensional pointer array of the same size as the `ndim` dimensional array.
    !>  </ol>
    !>  See the examples below.
    !>
    !>  \devnote
    !>  The logic behind allowing a pair of time series data `x` and `y` is to allow fast computation of the mean of the pair and the weight sums in one loop.<br>
    !>  This pattern occurs frequently in the computation of [correlation coefficients](@ref pm_sampleCor).<br>
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setVarMean](@ref pm_sampleVar::setVarMean)<br>
    !>
    !>  \example{setMean}
    !>  \include{lineno} example/pm_sampleMean/setMean/main.F90
    !>  \compilef{setMean}
    !>  \output{setMean}
    !>  \include{lineno} example/pm_sampleMean/setMean/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setMean_dim1_vs_dim2, The runtime performance of [setMean](@ref pm_sampleMean::setMean) along different sample dimensions.}
    !>  \include{lineno} benchmark/pm_sampleMean/setMean_dim1_vs_dim2/main.F90
    !>  \compilefb{setMean_dim1_vs_dim2}
    !>  \postprocb{setMean_dim1_vs_dim2}
    !>  \include{lineno} benchmark/pm_sampleMean/setMean_dim1_vs_dim2/main.py
    !>  \visb{setMean_dim1_vs_dim2}
    !>  \image html benchmark/pm_sampleMean/setMean_dim1_vs_dim2/benchmark.setMean_dim1_vs_dim2.runtime.png width=1000
    !>  \image html benchmark/pm_sampleMean/setMean_dim1_vs_dim2/benchmark.setMean_dim1_vs_dim2.runtime.ratio.png width=1000
    !>  \moralb{setMean_dim1_vs_dim2}
    !>      -#  The procedures under the generic interface [setMean](@ref pm_sampleMean::setMean) can compute the covariance under two different sample axes.<br>
    !>      -#  Recall that C is a row-major language while Fortran is a column-major language.<br>
    !>      -#  As such, one would expect the computations for a sample whose observations are stored along the second axis would be faster in the Fortran programming language.<br>
    !>      -#  However, such an expectation does not appear to hold at all times and appears to depend significantly on the computing architecture and the number of data attributes involved.<br>
    !>      -#  The higher the number of data attributes, the more likely the computations along the second axis of `sample` will be faster.<br>
    !>      -#  Note that for small number of data attributes, the computations along the second data axis involve a small
    !>          loop that has significant computational cost due to the implicit branching involved in the loop.<br>
    !>
    !>  \test
    !>  [test_pm_sampleMean](@ref test_pm_sampleMean)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! ALL XY

    interface setMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_CK5(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        complex(TKG)    , intent(out)                               :: mean(2)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_CK4(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        complex(TKG)    , intent(out)                               :: mean(2)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_CK3(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        complex(TKG)    , intent(out)                               :: mean(2)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_CK2(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        complex(TKG)    , intent(out)                               :: mean(2)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_CK1(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        complex(TKG)    , intent(out)                               :: mean(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_RK5(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(out)                               :: mean(2)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_RK4(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(out)                               :: mean(2)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_RK3(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(out)                               :: mean(2)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_RK2(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(out)                               :: mean(2)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WNO_XY_RK1(mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(out)                               :: mean(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_CK5(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_CK4(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_CK3(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_CK2(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_CK1(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_RK5(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_RK4(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_RK3(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_RK2(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WTI_XY_RK1(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_CK5(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_CK4(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_CK3(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_CK2(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_CK1(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_RK5(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_RK4(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_RK3(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_RK2(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WTR_XY_RK1(mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: x(:), y(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean(2)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ALL D1

    interface setMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_CK5(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_CK4(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_CK3(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_CK2(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_CK1(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_RK5(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_RK4(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_RK3(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_RK2(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WNO_D1_RK1(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_CK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_CK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_CK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_CK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_CK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_RK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_RK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_RK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_RK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WTI_D1_RK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_CK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_CK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_CK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_CK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_CK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_RK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_RK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_RK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_RK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WTR_D1_RK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ALL D2

    interface setMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_CK5(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_CK4(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_CK3(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_CK2(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_CK1(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_RK5(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_RK4(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_RK3(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_RK2(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WNO_D2_RK1(mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_CK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_CK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_CK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_CK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_CK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_RK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_RK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_RK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_RK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WTI_D2_RK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_CK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_CK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_CK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_CK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_CK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_RK5(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_RK4(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_RK3(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_RK2(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanALL_WTR_D2_RK1(mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanALL_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DIM D1

    interface setMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_CK5(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_CK4(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_CK3(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_CK2(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_CK1(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(out)                               :: mean
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_RK5(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_RK4(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_RK3(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_RK2(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanDIM_WNO_D1_RK1(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(out)                               :: mean
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_CK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_CK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_CK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_CK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_CK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_RK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_RK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_RK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_RK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanDIM_WTI_D1_RK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_CK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_CK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_CK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_CK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_CK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_RK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_RK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_RK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_RK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanDIM_WTR_D1_RK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DIM D2

    interface setMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_CK5(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_CK4(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_CK3(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_CK2(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_CK1(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_RK5(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_RK4(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_RK3(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_RK2(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanDIM_WNO_D2_RK1(mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_CK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_CK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_CK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_CK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_CK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_RK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_RK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_RK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_RK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanDIM_WTI_D2_RK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        integer(IK)     , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_CK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_CK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_CK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_CK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_CK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_RK5(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_RK4(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_RK3(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_RK2(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanDIM_WTR_D2_RK1(mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanDIM_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)                               :: weisum
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (weighted) merged mean of a sample resulting from the merger of two separate (weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampelMean](@ref pm_sampleMean) for more information and definition online updating of sample mean.<br>
    !>
    !>  \param[in]      meanB   :   The input object of the same type and kind and rank and shape as the input argument `meanA`,
    !>                              containing the mean of the second sample that must be merged with the first sample.<br>
    !>  \param[in]      meanA   :   The input scalar or `contiguous` vector of shape `(1:ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the mean of the first sample.<br>
    !>  \param[in]      fracA   :   The input scalar of type `real` of the same kind as kind of `meanA`,
    !>                              containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                              If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>
    !>  \return
    !>  `meanMerged`            :   The output object of the same type and kind and rank and shape as `meanA`,
    !>                              containing the mean of the sample resulting form the merger of the two samples.<br>
    !>
    !>  \interface{getMeanMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleMean, only: getMeanMerged
    !>
    !>      ! univariate sample.
    !>
    !>      meanMerged = getMeanMerged(meanB, meanA, fracA)
    !>
    !>      ! ndim-dimensional sample.
    !>
    !>      meanMerged(1:ndim) = getMeanMerged(meanB(1:ndim), meanA(1:ndim), fracA)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(meanB) == size(meanA)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>
    !>  \example{setMeanMerged}
    !>  \include{lineno} example/pm_sampleMean/setMeanMerged/main.F90
    !>  \compilef{setMeanMerged}
    !>  \output{setMeanMerged}
    !>  \include{lineno} example/pm_sampleMean/setMeanMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleMean](@ref test_pm_sampleMean)
    !>
    !>  \final{getMeanMerged}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! D0

    interface getMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanMergedNew_D0_CK5(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)                                                :: meanMerged
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanMergedNew_D0_CK4(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)                                                :: meanMerged
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanMergedNew_D0_CK3(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)                                                :: meanMerged
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanMergedNew_D0_CK2(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)                                                :: meanMerged
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanMergedNew_D0_CK1(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)                                                :: meanMerged
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanMergedNew_D0_RK5(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)                                                   :: meanMerged
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanMergedNew_D0_RK4(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)                                                   :: meanMerged
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanMergedNew_D0_RK3(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)                                                   :: meanMerged
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanMergedNew_D0_RK2(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)                                                   :: meanMerged
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanMergedNew_D0_RK1(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)                                                   :: meanMerged
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1

    interface getMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMeanMergedNew_D1_RK5(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)                                                   :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMeanMergedNew_D1_RK4(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)                                                   :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMeanMergedNew_D1_RK3(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)                                                   :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMeanMergedNew_D1_RK2(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)                                                   :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMeanMergedNew_D1_RK1(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)                                                   :: meanMerged(size(meanA, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMeanMergedNew_D1_CK5(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)                                                :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMeanMergedNew_D1_CK4(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)                                                :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMeanMergedNew_D1_CK3(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)                                                :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMeanMergedNew_D1_CK2(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)                                                :: meanMerged(size(meanA, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMeanMergedNew_D1_CK1(meanB, meanA, fracA) result(meanMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMergedNew_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)                                                :: meanMerged(size(meanA, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the (weighted) merged mean of a sample resulting from the merger of two separate (weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampelMean](@ref pm_sampleMean) for more information and definition online updating of sample mean.<br>
    !>
    !>  \param[out]     meanMerged  :   The output object of the same type and kind and rank and shape as `meanA`,
    !>                                  containing the mean of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged mean will be written to the argument `meanB`.)
    !>  \param[inout]   meanB       :   The input or input/output object of the same type and kind and rank and shape as the input argument `meanA`,
    !>                                  containing the mean of the second sample that must be merged with the first sample.<br>
    !>                                  If the input argument `meanMerged` is missing, then `meanB` contains the merged mean of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `meanB` remain intact upon return.<br>
    !>  \param[in]      meanA       :   The input scalar or `contiguous` vector of shape `(1:ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the mean of the first sample.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as kind of `meanA`,
    !>                                  containing the sum of the weights of all points in sample \f$A\f$ divided by the sum of the weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>
    !>  \interface{setMeanMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleMean, only: setMeanMerged
    !>
    !>      ! univariate sample.
    !>
    !>      call setMeanMerged(            meanB, meanA, fracA)
    !>      call setMeanMerged(meanMerged, meanB, meanA, fracA)
    !>
    !>      ! ndim-dimensional sample.
    !>
    !>      call setMeanMerged(                    meanB(1:ndim), meanA(1:ndim), fracA)
    !>      call setMeanMerged(meanMerged(1:ndim), meanB(1:ndim), meanA(1:ndim), fracA)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(meanB) == size(meanA)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(meanMerged) == size(meanA)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>
    !>  \example{setMeanMerged}
    !>  \include{lineno} example/pm_sampleMean/setMeanMerged/main.F90
    !>  \compilef{setMeanMerged}
    !>  \output{setMeanMerged}
    !>  \include{lineno} example/pm_sampleMean/setMeanMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleMean](@ref test_pm_sampleMean)
    !>
    !>  \final{setMeanMerged}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! New_D0

    interface setMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanMergedNew_D0_CK5(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanMergedNew_D0_CK4(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanMergedNew_D0_CK3(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanMergedNew_D0_CK2(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanMergedNew_D0_CK1(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanMergedNew_D0_RK5(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanMergedNew_D0_RK4(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanMergedNew_D0_RK3(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanMergedNew_D0_RK2(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanMergedNew_D0_RK1(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! New_D1

    interface setMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanMergedNew_D1_CK5(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanMergedNew_D1_CK4(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanMergedNew_D1_CK3(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanMergedNew_D1_CK2(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanMergedNew_D1_CK1(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:), meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanMergedNew_D1_RK5(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanMergedNew_D1_RK4(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanMergedNew_D1_RK3(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanMergedNew_D1_RK2(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanMergedNew_D1_RK1(meanMerged, meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedNew_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:), meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_D0

    interface setMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanMergedOld_D0_CK5(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(inout)                             :: meanB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanMergedOld_D0_CK4(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(inout)                             :: meanB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanMergedOld_D0_CK3(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(inout)                             :: meanB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanMergedOld_D0_CK2(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(inout)                             :: meanB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanMergedOld_D0_CK1(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(inout)                             :: meanB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanMergedOld_D0_RK5(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(inout)                             :: meanB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanMergedOld_D0_RK4(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(inout)                             :: meanB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanMergedOld_D0_RK3(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(inout)                             :: meanB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanMergedOld_D0_RK2(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(inout)                             :: meanB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanMergedOld_D0_RK1(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(inout)                             :: meanB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_D1

    interface setMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMeanMergedOld_D1_CK5(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMeanMergedOld_D1_CK4(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMeanMergedOld_D1_CK3(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMeanMergedOld_D1_CK2(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMeanMergedOld_D1_CK1(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMeanMergedOld_D1_RK5(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMeanMergedOld_D1_RK4(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMeanMergedOld_D1_RK3(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMeanMergedOld_D1_RK2(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMeanMergedOld_D1_RK1(meanB, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMeanMergedOld_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(inout) , contiguous                :: meanB(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleMean ! LCOV_EXCL_LINE
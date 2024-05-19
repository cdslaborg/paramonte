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
!>  This module contains classes and procedures for computing properties related to the **auto correlation time (ACT)** of random sequences.
!>
!>  \details
!>  The integrated autocorrelation time (ACT or IAT or IAC) frequently appears in the calculation of
!>  the effective sample size of the outputs of [Markov Chain Monte Carlo samplers](@ref pm_sampling).<br>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>  [pm_sampleACT](@ref pm_sampleACT)<br>
!>  [pm_sampleCCF](@ref pm_sampleCCF)<br>
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>  [pm_sampleCov](@ref pm_sampleCov)<br>
!>  [pm_sampleVar](@ref pm_sampleVar)<br>
!>  [pm_sampleConv](@ref pm_sampleConv)<br>
!>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
!>  [pm_sampleMean](@ref pm_sampleMean)<br>
!>  [pm_sampleNorm](@ref pm_sampleNorm)<br>
!>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
!>  [pm_sampleScale](@ref pm_sampleScale)<br>
!>  [pm_sampleShift](@ref pm_sampleShift)<br>
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_sampleAffinity](@ref pm_sampleAffinity)<br>
!>
!>  \test
!>  [test_pm_sampleACT](@ref test_pm_sampleACT)
!>
!>  \final{pm_sampleACT}
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 01:45 AM, August 21, 2018, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleACT

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleACT"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, abstract :: act_type
    end type

    type, extends(act_type) :: cumSum_type
        real :: signif = 0. !<  \public The autocorrelation significance in units of standard deviation below which the autocorrelation is considered noise.<br>
    end type

    type, extends(act_type) :: cumSumMax_type
    end type

    type, extends(act_type) :: batchMeans_type
        integer(IK) :: size = 0_IK
    end type

    type, extends(act_type) :: batchMeansMax_type
        integer(IK) :: smin = 2_IK          !<  \public The scalar `integer` of default kind \IK, containing the minimum batch size to be considered.
        integer(IK) :: smax = huge(0_IK)    !<  \public The scalar `integer` of default kind \IK, containing the maximum batch size to be considered.
        integer(IK) :: step = 1_IK          !<  \public The scalar `integer` of default kind \IK, containing the jump in batch size to be considered.
    end type

    type(cumSum_type), parameter :: cumSum = cumSum_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cumSum
#endif
    type(cumSumMax_type), parameter :: cumSumMax = cumSumMax_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cumSumMax
#endif
    type(batchMeans_type), parameter :: batchMeans = batchMeans_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: batchMeans
#endif
    type(batchMeansMax_type ), parameter :: batchMeansMax = batchMeansMax_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: batchMeansMax
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the integrated auto-correlation time (ACT) of the discrete signal.<br>
    !>
    !>  \details
    !>  This generic interface a flexible wrapper for the lower-level potentially faster generic interface [setACF](@ref pm_sampleCCF::setACF).<br>
    !>  See the documentation of the parent module [pm_sampleACT](@ref pm_sampleACT) for algorithmic details and auto-correlation definition.<br>
    !>
    !>  \param[in]  seq         :   The input `contiguous` vector of arbitrary size (of minimum `2`) of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sequence auto-correlation time must be computed.<br>
    !>  \param[in]  weight      :   The input vector of type `integer` of default kind \IK of the same size as the input `seq`,
    !>                              containing the weights of the corresponding elements of the input `seq`.<br>
    !>                              (**optional**, default = `[(1, i = 1, size(seq))]`)
    !>  \param[in]  method      :   The input scalar constant that can be:<br>
    !>                              <ol>
    !>                                  <li>    the constant [batchMeans](@ref pm_sampleACT::batchMeans) or an object of type [batchMeans_type(size = size)](@ref pm_sampleACT::batchMeans_type),
    !>                                          implying that the returned ACT is computed using the BatchMeans method with the optionally specified `size`.<br>
    !>                                          If the `size` component of the input object of type [batchMeans_type(size = size)](@ref pm_sampleACT::batchMeans_type) is unset or `0`,
    !>                                          an appropriate default based on the length of the input sequence will be determined and used.<br>
    !>                                          This option offers a reasonable tradeoff between good estimate and computational efficiency.<br>
    !>                                          However, the results depend heavily on the specified batch size, hence, still underestimate the actual ACT.<br>
    !>                                  <li>    the constant [batchMeansMax](@ref pm_sampleACT::batchMeansMax) or an object of type [batchMeans_type(smin, smax, step)](@ref pm_sampleACT::batchMeans_type),
    !>                                          implying that the returned ACT is computed as the maximum of the BatchMeans method for the specified range of
    !>                                          batch sizes `getRange(batchMeans_type%%smin, batchMeans_type%%smax, batchMeans_type%%step).<br>
    !>                                          If the components of the input object of type [batchMeans_type(smin, smax, step)](@ref pm_sampleACT::batchMeans_type) are unset (default),
    !>                                          then `smin` and `step` will be both set to `1` and `smax` is set to `size(seq) / 2`.<br>
    !>                                          This option will likely yield the least biased results but will be computationally significantly more demanding.<br>
    !>                                  <li>    the constant [cumSum](@ref pm_sampleACT::cumSum) or an object of type [cumSum_type(signif)](@ref pm_sampleACT::cumSum_type),
    !>                                          implying that the returned ACT is simply the cumulative sum of auto-correlations of the input sequence at
    !>                                          all possible lags: `act = 1 + 2 * sum(getACT(seq))`.<br>
    !>                                          This option is merely provided for exploratory and education reasons.<br>
    !>                                  <li>    the constant [cumSumMax](@ref pm_sampleACT::cumSumMax) or an object of type [cumSumMax_type(smin, smax, step)](@ref pm_sampleACT::cumSumMax_type),
    !>                                          implying that the returned ACT is computed as the maximum of the cumulative sum of auto-correlations of the input sequence at all possible lags.<br>
    !>                                          This option is merely provided for exploratory and education reasons.<br>
    !>                                          This option is merely provided for exploratory and education reasons.<br>
    !>                              </ol>
    !>                              (**optional**, default = [batchMeans_type(size = size(seq)**0.66)](@ref pm_sampleACT::batchMeans).)
    !>
    !>  \return
    !>  `act`                   :   The output scalar of the same type and kind as that of the input sequence,
    !>                              containing the estimated integrated autocorrelation time of the input sequence.<br>
    !>
    !>  \interface{getACT}
    !>  \code{.F90}
    !>
    !>      use pm_sampleACT, only: getACT
    !>
    !>      act = getACT(seq(:), method = method)
    !>      act = getACT(seq(:), weight, method = method)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < method%smin` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < method%signif` must hold for the corresponding input arguments.<br>
    !>  The condition `size(weight) == size(seq)` must hold for the corresponding input arguments.<br>
    !>  The condition `method%smin <= method%smax <= size(seq, 1) / 2` must hold for the corresponding input arguments.<br>
    !>  The condition `method%size * 2 < sum(weight)` must hold for the corresponding input arguments, otherwise `sum(weight)` will be returned as the `act`.<br>
    !>  The condition `1 < size(seq)` must hold for the corresponding input arguments, otherwise `act` will be set to `lenseq` on return.<br>
    !>  The condition `0 < method%step` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < method%size` must hold for the corresponding input arguments.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getACT](@ref pm_sampleACT::getACT)<br>
    !>  [getACF](@ref pm_sampleCCF::getACF)<br>
    !>  [setACF](@ref pm_sampleCCF::setACF)<br>
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{getACT}
    !>  \include{lineno} example/pm_sampleACT/getACT/main.F90
    !>  \compilef{getACT}
    !>  \output{getACT}
    !>  \include{lineno} example/pm_sampleACT/getACT/main.out.F90
    !>  \postproc{getACT}
    !>  \include{lineno} example/pm_sampleACT/getACT/main.py
    !>  \vis{getACT}
    !>  \image html pm_sampleACT/getACT/getACT.duluth.thinned.png width=700
    !>  \image html pm_sampleACT/getACT/getACT.duluth.batchMeans.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleACT](@ref test_pm_sampleACT)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      getACT_CSD_WTI_D1_CK5()
    !>         ||| ||| ||| || |||
    !>         ||| ||| ||| || The type and kind parameters of the input sequence.
    !>         ||| ||| ||| The dimension of the input sequence `seq`.
    !>         ||| ||| Weight type: WTTI => integer weight, ONE => no weight.
    !>         ||| The method: BMD => BatchMeans Default, BMM => BatchMeans Maximum, CSD => CumSum Default, CSM => CumSum Maximum
    !>         ACF: Cross-Correlation Function.
    !>  \endcode
    !>
    !>  \todo
    !>  \phigh
    !>  The performance of the implementations of this generic interface can be improved by minimizing allocations and converting function calls to subroutine calls.<br>
    !>
    !>  \final{getACT}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! D1

    interface getACT

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_DEF_ONE_D1_RK5(seq) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_ONE_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_DEF_ONE_D1_RK4(seq) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_ONE_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_DEF_ONE_D1_RK3(seq) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_ONE_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_DEF_ONE_D1_RK2(seq) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_ONE_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_DEF_ONE_D1_RK1(seq) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_ONE_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_DEF_WTI_D1_RK5(seq, weight) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_WTI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_DEF_WTI_D1_RK4(seq, weight) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_WTI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_DEF_WTI_D1_RK3(seq, weight) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_WTI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_DEF_WTI_D1_RK2(seq, weight) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_WTI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_DEF_WTI_D1_RK1(seq, weight) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_DEF_WTI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_CSD_ONE_D1_RK5(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_ONE_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(cumSum_type)       , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_CSD_ONE_D1_RK4(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_ONE_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(cumSum_type)       , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_CSD_ONE_D1_RK3(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_ONE_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(cumSum_type)       , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_CSD_ONE_D1_RK2(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_ONE_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(cumSum_type)       , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_CSD_ONE_D1_RK1(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_ONE_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(cumSum_type)       , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_CSD_WTI_D1_RK5(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_WTI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(cumSum_type)       , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_CSD_WTI_D1_RK4(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_WTI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(cumSum_type)       , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_CSD_WTI_D1_RK3(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_WTI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(cumSum_type)       , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_CSD_WTI_D1_RK2(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_WTI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(cumSum_type)       , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_CSD_WTI_D1_RK1(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSD_WTI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(cumSum_type)       , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_CSM_ONE_D1_RK5(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_ONE_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(cumSumMax_type)    , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_CSM_ONE_D1_RK4(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_ONE_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(cumSumMax_type)    , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_CSM_ONE_D1_RK3(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_ONE_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(cumSumMax_type)    , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_CSM_ONE_D1_RK2(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_ONE_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(cumSumMax_type)    , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_CSM_ONE_D1_RK1(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_ONE_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(cumSumMax_type)    , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_CSM_WTI_D1_RK5(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_WTI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(cumSumMax_type)    , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_CSM_WTI_D1_RK4(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_WTI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(cumSumMax_type)    , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_CSM_WTI_D1_RK3(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_WTI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(cumSumMax_type)    , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_CSM_WTI_D1_RK2(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_WTI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(cumSumMax_type)    , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_CSM_WTI_D1_RK1(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_CSM_WTI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(cumSumMax_type)    , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_BMD_ONE_D1_RK5(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_ONE_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(batchMeans_type)   , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_BMD_ONE_D1_RK4(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_ONE_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(batchMeans_type)   , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_BMD_ONE_D1_RK3(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_ONE_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(batchMeans_type)   , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_BMD_ONE_D1_RK2(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_ONE_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(batchMeans_type)   , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_BMD_ONE_D1_RK1(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_ONE_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(batchMeans_type)   , intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_BMD_WTI_D1_RK5(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_WTI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(batchMeans_type)   , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_BMD_WTI_D1_RK4(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_WTI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(batchMeans_type)   , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_BMD_WTI_D1_RK3(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_WTI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(batchMeans_type)   , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_BMD_WTI_D1_RK2(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_WTI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(batchMeans_type)   , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_BMD_WTI_D1_RK1(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMD_WTI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(batchMeans_type)   , intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_BMM_ONE_D1_RK5(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_ONE_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(batchMeansMax_type), intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_BMM_ONE_D1_RK4(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_ONE_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(batchMeansMax_type), intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_BMM_ONE_D1_RK3(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_ONE_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(batchMeansMax_type), intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_BMM_ONE_D1_RK2(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_ONE_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(batchMeansMax_type), intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_BMM_ONE_D1_RK1(seq, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_ONE_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(batchMeansMax_type), intent(in)                            :: method
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACT_BMM_WTI_D1_RK5(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_WTI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(batchMeansMax_type), intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK4_ENABLED
    module function getACT_BMM_WTI_D1_RK4(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_WTI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(batchMeansMax_type), intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK3_ENABLED
    module function getACT_BMM_WTI_D1_RK3(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_WTI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(batchMeansMax_type), intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK2_ENABLED
    module function getACT_BMM_WTI_D1_RK2(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_WTI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(batchMeansMax_type), intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

#if RK1_ENABLED
    module function getACT_BMM_WTI_D1_RK1(seq, weight, method) result(act)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACT_BMM_WTI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(batchMeansMax_type), intent(in)                            :: method
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(RKG)               , intent(in)    , contiguous            :: seq(:)
        real(RKG)                                                       :: act
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getACT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleACT ! LCOV_EXCL_LINE
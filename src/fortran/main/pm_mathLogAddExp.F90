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
!>  This module contains procedures and generic interfaces
!>  for adding two real or complex values without causing overflow or underflow.
!>
!>  \see
!>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
!>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
!>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
!>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
!>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
!>
!>  \finmain
!>
!>  \test
!>  [test_pm_mathLogAddExp](@ref test_pm_mathLogAddExp)
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathLogAddExp

    use pm_kind, only: SK, IK, LK
    use pm_control, only: selection_type, selection

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathLogAddExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the logarithm of the sum of the exponential of two input logarithmic values **robustly (without causing overflow)**.
    !>
    !>  \param[in]  smaller     :   The input scalar, or array of the same rank, shape, and size as other array-like input arguments, of either <br>
    !>                              <ul>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              representing the smaller of the two logarithmic values to be added.<br>
    !>                              If it is `complex`, then the `real` component of the `complex` value
    !>                              must be smaller than the `real` component of the (`larger`) `complex` value.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> the input argument `larger` is also present.)
    !>  \param[in]  larger      :   The input scalar, or array of the same rank, shape, and size as other array-like input arguments, of the same
    !>                              type and kind as the input argument `smaller`, representing the larger of the two logarithmic values to be added.<br>
    !>                              If it is `complex`, then the `real` component of the `complex` value
    !>                              must be larger than the `real` component of the (`larger`) `complex` value.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> the input argument `smaller` is also present.)
    !>  \param[in]  minMax      :   The input vector of length `2` of either <br>
    !>                              <ul>
    !>                                  <li>    type `complex` of kind \CKALL or, <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              representing the smaller and the larger of the two logarithmic values to be added as `minMax(1)` and `minMax(2)`.<br>
    !>                              If it is `complex`, then the `real` component of the `complex` value must be
    !>                              smaller than the `real` component of the (`larger`) `complex` value.<br>
    !>                              If the `smaller` and `larger` values are not known a priori, the vector `minMax`
    !>                              can be instead constructed by calling [getMinMax](@ref pm_mathMinMax::getMinMax).<br>
    !>                              (**optional**, it must be present <b>if and only if</b> the input arguments `smaller` and `larger` are missing.)
    !>  \param      control     :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [selection](@ref pm_control::selection) or equivalently,
    !>                                          an object of type [selection_type](@ref pm_control::selection_type).<br>
    !>                                          Specifying this option enables the runtime checks for underflow occurrence via branching and dynamic dispatch.<br>
    !>                                          Enabling this option can aid runtime efficiency when the `smaller` to `larger` ratio is expected to cause underflow.<br>
    !>                                          In such cases, the exponentiation is avoided if `control =` [selection](@ref pm_control::selection), leading to faster
    !>                                          runtime by avoiding exponentiation since it is highly expensive (on the order of ~200 CPU cycles).<br>
    !>                                          See [the relevant benchmark here](#benchmark-getLogAddExp).<br>
    !>                              </ol>
    !>                              The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                              (**optional**. If missing, a sequence control flow will be assumed without dynamic dispatch.)
    !>
    !>  \return
    !>  `logAddExp`             :   The output scalar of the same type and kind as the input `larger` and `smaller` representing
    !>                              the logarithm of the sum of the exponential of two input logarithmic values computed robustly.
    !>
    !>  \interface{getLogAddExp}
    !>  \code{.F90}
    !>
    !>      use pm_mathLogAddExp, only: getLogAddExp, selection
    !>
    !>      logAddExp = getLogAddExp(smaller, larger)
    !>      logAddExp = getLogAddExp(smaller, larger, control)
    !>
    !>      logAddExp = getLogAddExp(minMax(1:2))
    !>      logAddExp = getLogAddExp(minMax(1:2), control)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  As the argument names suggest, the value of the input argument `larger` must be larger than or equal to the value of the input argument `smaller`.<br>
    !>  When the input arguments are of type `complex`, this condition must hold for the `real` components of the numbers.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when the input arguments `smaller` and `larger` are present.<br>
    !>
    !>  \note
    !>  The smaller and larger input values can be quickly obtained by a one-line call to [setMinMax](@ref pm_mathMinMax::setMinMax).<br>
    !>  Alternatively, the smaller and larger values can be sorted in a `minMax(1:2)` array by calling [getMinMax](@ref pm_mathMinMax::getMinMax)
    !>  and passing the result directly to [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp).<br>
    !>
    !>  \see
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
    !>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
    !>
    !>  \example{getLogAddExp}
    !>  \include{lineno} example/pm_mathLogAddExp/getLogAddExp/main.F90
    !>  \compilef{getLogAddExp}
    !>  \output{getLogAddExp}
    !>  \include{lineno} example/pm_mathLogAddExp/getLogAddExp/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getLogAddExp, The effects of `control` argument on runtime efficiency}
    !>      The following program compares the runtime performance of [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)
    !>      algorithm with and without checking for underflows.
    !>  \include{lineno} benchmark/pm_mathLogAddExp/getLogAddExp/main.F90
    !>  \compilefb{getLogAddExp}
    !>  \postprocb{getLogAddExp}
    !>  \include{lineno} benchmark/pm_mathLogAddExp/getLogAddExp/main.py
    !>  \visb{getLogAddExp}
    !>  \image html benchmark/pm_mathLogAddExp/getLogAddExp/benchmark.getLogAddExp.normal.png width=1000
    !>  \image html benchmark/pm_mathLogAddExp/getLogAddExp/benchmark.getLogAddExp.underflow.png width=1000
    !>  \moralb{getLogAddExp}
    !>      -#  If the ratio of the `smaller` to `larger` input arguments causes frequent underflows,
    !>          then setting `control =` [selection](@ref pm_control::selection) to allow branching will likely result in a faster runtime.<br>
    !>          Conversely, if the divisions are **not** expected to cause any or too many underflows,
    !>          then removing the input argument `control` when calling [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp) will lead
    !>          to a better runtime performance (at the expense of occasional expensive but redundant exponentiation operations).<br>
    !>      -#  If the procedure [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp) is to be
    !>          called billions of times, then it would make sense to manually inline the procedure implementation in your code as
    !>          procedure call will have a non-negligible performance overhead.<br>
    !>
    !>  \test
    !>  [test_pm_mathLogAddExp](@ref test_pm_mathLogAddExp)
    !>
    !>  \finmain{getLogAddExp}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface getLogAddExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getLogAddExpSeqSL_CK5(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getLogAddExpSeqSL_CK4(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getLogAddExpSeqSL_CK3(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getLogAddExpSeqSL_CK2(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getLogAddExpSeqSL_CK1(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogAddExpSeqSL_RK5(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogAddExpSeqSL_RK4(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogAddExpSeqSL_RK3(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogAddExpSeqSL_RK2(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogAddExpSeqSL_RK1(smaller, larger) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqSL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getLogAddExpSelSL_CK5(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getLogAddExpSelSL_CK4(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getLogAddExpSelSL_CK3(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getLogAddExpSelSL_CK2(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getLogAddExpSelSL_CK1(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogAddExpSelSL_RK5(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogAddExpSelSL_RK4(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogAddExpSelSL_RK3(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogAddExpSelSL_RK2(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogAddExpSelSL_RK1(smaller, larger, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelSL_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLogAddExpSeqMM_CK5(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogAddExpSeqMM_CK4(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogAddExpSeqMM_CK3(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogAddExpSeqMM_CK2(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogAddExpSeqMM_CK1(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogAddExpSeqMM_RK5(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogAddExpSeqMM_RK4(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogAddExpSeqMM_RK3(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogAddExpSeqMM_RK2(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogAddExpSeqMM_RK1(minMax) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSeqMM_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLogAddExpSelMM_CK5(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogAddExpSelMM_CK4(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogAddExpSelMM_CK3(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogAddExpSelMM_CK2(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogAddExpSelMM_CK1(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logAddExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogAddExpSelMM_RK5(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogAddExpSelMM_RK4(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogAddExpSelMM_RK3(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogAddExpSelMM_RK2(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogAddExpSelMM_RK1(minMax, control) result(logAddExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogAddExpSelMM_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logAddExp
    end function
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

end module pm_mathLogAddExp ! LCOV_EXCL_LINE
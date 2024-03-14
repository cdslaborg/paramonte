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
!>  This module contains procedures and generic interfaces for subtracting two real or complex values without causing overflow or underflow.
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
!>  [test_pm_mathLogSubExp](@ref test_pm_mathLogSubExp)
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathLogSubExp

    use pm_kind, only: SK, IK, LK
    use pm_control, only: selection_type, selection

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathLogSubExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the subtraction of the exponential of the smaller
    !>  input logarithmic value from the larger input logarithmic value **robustly (without causing overflow)**.
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
    !>                              If it is `complex`, then the `real` component of the `complex` value must be smaller
    !>                              than the `real` component of the (`larger`) `complex` value.<br>
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
    !>                                          See [the relevant benchmark here](#benchmark-getLogSubExp).<br>
    !>                              </ol>
    !>                              The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                              (**optional**. If missing, a sequence control flow will be assumed without dynamic dispatch.)
    !>
    !>  \return
    !>  `logSubExp`             :   The output scalar of the same type and kind as the input `larger` and `smaller`
    !>                              containing the natural logarithm of the subtraction of the exponential of the larger
    !>                              input logarithmic value from the smaller input logarithmic value, computed robustly.<br>
    !>
    !>  \interface{getLogSubExp}
    !>  \code{.F90}
    !>
    !>      use pm_mathLogSubExp, only: getLogSubExp
    !>
    !>      logSubExp = getLogSubExp(smaller, larger)
    !>      logSubExp = getLogSubExp(smaller, larger, control)
    !>      logSubExp = getLogSubExp(minMax(1:2))
    !>      logSubExp = getLogSubExp(minMax(1:2), control)
    !>      !
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
    !>  and passing the result directly to [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp).<br>
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
    !>  \example{getLogSubExp}
    !>  \include{lineno} example/pm_mathLogSubExp/getLogSubExp/main.F90
    !>  \compilef{getLogSubExp}
    !>  \output{getLogSubExp}
    !>  \include{lineno} example/pm_mathLogSubExp/getLogSubExp/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getLogSubExp, The effects of the `control` argument on runtime efficiency}
    !>      The following program compares the runtime performance of [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)
    !>      algorithm with and without checking for underflows.<br>
    !>  \include{lineno} benchmark/pm_mathLogSubExp/getLogSubExp/main.F90
    !>  \compilefb{getLogSubExp}
    !>  \postprocb{getLogSubExp}
    !>  \include{lineno} benchmark/pm_mathLogSubExp/getLogSubExp/main.py
    !>  \visb{getLogSubExp}
    !>  \image html benchmark/pm_mathLogSubExp/getLogSubExp/benchmark.getLogSubExp.normal.png width=1000
    !>  \image html benchmark/pm_mathLogSubExp/getLogSubExp/benchmark.getLogSubExp.underflow.png width=1000
    !>  \moralb{getLogSubExp}
    !>      -#  If the ratio of the `smaller` to `larger` input arguments causes frequent underflows,
    !>          then setting `control =` [selection](@ref pm_control::selection) to allow branching will likely result in a faster runtime.<br>
    !>          Conversely, if the divisions are **not** expected to cause any or too many underflows, 
    !>          then removing the input argument `control` when calling [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp) will lead
    !>          to a better runtime performance (at the expense of occasional expensive but redundant exponentiation operations).<br>
    !>      -#  If the procedure [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp) is to be
    !>          called billions of times, then it would make sense to manually inline the procedure implementation in your code as
    !>          procedure call and processing of optional arguments will have a non-negligible performance overhead.<br>
    !>
    !>  \test
    !>  [test_pm_mathLogSubExp](@ref test_pm_mathLogSubExp)
    !>
    !>  \finmain{getLogSubExp}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface getLogSubExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getLogSubExpSeqSL_CK5(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getLogSubExpSeqSL_CK4(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getLogSubExpSeqSL_CK3(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getLogSubExpSeqSL_CK2(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getLogSubExpSeqSL_CK1(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogSubExpSeqSL_RK5(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogSubExpSeqSL_RK4(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogSubExpSeqSL_RK3(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogSubExpSeqSL_RK2(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogSubExpSeqSL_RK1(smaller, larger) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqSL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getLogSubExpSelSL_CK5(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getLogSubExpSelSL_CK4(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getLogSubExpSelSL_CK3(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getLogSubExpSelSL_CK2(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getLogSubExpSelSL_CK1(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: smaller, larger
        complex(CKC)                        :: logSubExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogSubExpSelSL_RK5(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogSubExpSelSL_RK4(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogSubExpSelSL_RK3(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogSubExpSelSL_RK2(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogSubExpSelSL_RK1(smaller, larger, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelSL_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: smaller, larger
        real(RKC)                           :: logSubExp
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
    PURE module function getLogSubExpSeqMM_CK5(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogSubExpSeqMM_CK4(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogSubExpSeqMM_CK3(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogSubExpSeqMM_CK2(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogSubExpSeqMM_CK1(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogSubExpSeqMM_RK5(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogSubExpSeqMM_RK4(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogSubExpSeqMM_RK3(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogSubExpSeqMM_RK2(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogSubExpSeqMM_RK1(minMax) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSeqMM_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLogSubExpSelMM_CK5(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogSubExpSelMM_CK4(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogSubExpSelMM_CK3(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogSubExpSelMM_CK2(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogSubExpSelMM_CK1(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(selection_type), intent(in)    :: control
        complex(CKC)        , intent(in)    :: minMax(2)
        complex(CKC)                        :: logSubExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogSubExpSelMM_RK5(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogSubExpSelMM_RK4(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogSubExpSelMM_RK3(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogSubExpSelMM_RK2(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogSubExpSelMM_RK1(minMax, control) result(logSubExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExpSelMM_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(selection_type), intent(in)    :: control
        real(RKC)           , intent(in)    :: minMax(2)
        real(RKC)                           :: logSubExp
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

end module pm_mathLogSubExp ! LCOV_EXCL_LINE
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
!>  This module contains the procedures and interfaces for computing the natural logarithm of the sum of exponentials the elements of an array.
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
!>  \test
!>  [test_pm_mathLogSumExp](@ref test_pm_mathLogSumExp)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathLogSumExp

    use pm_kind, only: SK, IK, LK
    use pm_control, only: selection_type, selection

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathLogSumExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the sum of the exponential of the input `array` robustly (without numerical overflow).
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either,
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL or, <br>
    !>                                  <li>    type `complex` of kind \CKALL, <br>
    !>                              </ul>
    !>                              whose natural logarithm of sum of the exponential of its elements will be generated and returned.<br>
    !>  \param[in]  maxArray    :   The input scalar of the same type and kind as the input `array` representing the maximum value in `array` (i.e., `maxArray = maxval(array)`).<br>
    !>                              If the input `array` is of `complex` type, then only `real` component must be considered in the computation of the maximum value.<br>
    !>                              Specifying this argument, if known a priori, will significantly aid the runtime performance.<br>
    !>                              (**optional**. It must be present if the input argument `control` is also present.)
    !>  \param      control     :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [selection](@ref pm_control::selection) or equivalently,
    !>                                          an object of type [selection_type](@ref pm_control::selection_type).<br>
    !>                                          Specifying this option enables the runtime checks for underflow occurrence via branching and dynamic dispatch.<br>
    !>                                          Enabling this option can aid runtime efficiency when the division of a significant number
    !>                                          of elements of `array` (for example, half or more) by `maxArray` causes underflow.<br>
    !>                                          This option avoids the computation of an exponentiation term, leading to better runtime efficiency.<br>
    !>                                          Note that exponentiation is highly expensive (on the order of ~200 CPU cycles).<br>
    !>                                          See [the relevant benchmark here](#benchmark-getLogSumExp).<br>
    !>                              </ol>
    !>                              The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                              (**optional**. If missing, a [sequence](@ref pm_control::sequence) control flow will be assumed without dynamic dispatch.)
    !>
    !>  \return
    !>  `logSumExp`            :    The output scalar of the same type and kind as the input `array` containing the natural
    !>                              logarithm of the sum of the exponential of the input `array` robustly (without numerical overflow).
    !>
    !>  \interface{getLogSumExp}
    !>  \code{.F90}
    !>
    !>      use pm_mathLogSumExp, only: getLogSumExp, selection
    !>
    !>      logSumExp = getLogSumExp(array(1:np), maxArray)
    !>      logSumExp = getLogSumExp(array(1:np), maxArray, control)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `real(maxArray) == maxval(real(array))` must hold.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
    !>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
    !>  [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp)<br>
    !>  [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)<br>
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [setCumSum](@ref pm_mathCumSum::setCumSum)<br>
    !>
    !>  \example{getLogSumExp}
    !>  \include{lineno} example/pm_mathLogSumExp/getLogSumExp/main.F90
    !>  \compilef{getLogSumExp}
    !>  \output{getLogSumExp}
    !>  \include{lineno} example/pm_mathLogSumExp/getLogSumExp/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getLogSumExp, The effects of `control` on runtime efficiency}
    !>      The following program compares the runtime performance of [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)
    !>      algorithm with and without checking for underflows.
    !>  \include{lineno} benchmark/pm_mathLogSumExp/getLogSumExp/main.F90
    !>  \compilefb{getLogSumExp}
    !>  \postprocb{getLogSumExp}
    !>  \include{lineno} benchmark/pm_mathLogSumExp/getLogSumExp/main.py
    !>  \visb{getLogSumExp}
    !>  \image html benchmark/pm_mathLogSumExp/getLogSumExp/benchmark.getLogSumExp.normal.png width=1000
    !>  \image html benchmark/pm_mathLogSumExp/getLogSumExp/benchmark.getLogSumExp.underflow.png width=1000
    !>  \moralb{getLogSumExp}
    !>      -#  If the input array has many (half the size of array or more) elements whose division by the `maxval(array)` causes underflow,
    !>          then setting `control =` [selection](@ref pm_control::selection) to allow branching will likely result in a faster runtime.<br>
    !>          Conversely, if the divisions are **not** expected to cause any or too many underflows,
    !>          then removing the input argument `control` when calling [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp) will lead
    !>          to a better runtime performance (at the expense of occasional expensive but redundant exponentiation operations).<br>
    !>      -#  Note the performance of manual computation against [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)
    !>          with or without setting the optional argument `control`.<br>
    !>
    !>  \test
    !>  [test_pm_mathLogSumExp](@ref test_pm_mathLogSumExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be expanded to include input arrays with `Weight`s.
    !>
    !>  \final{getLogSumExp}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface getLogSumExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLogSumExpDefSeq_CK5(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogSumExpDefSeq_CK4(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogSumExpDefSeq_CK3(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogSumExpDefSeq_CK2(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogSumExpDefSeq_CK1(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG)                            :: logSumExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogSumExpDefSeq_RK5(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogSumExpDefSeq_RK4(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogSumExpDefSeq_RK3(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogSumExpDefSeq_RK2(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogSumExpDefSeq_RK1(array) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpDefSeq_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)                               :: logSumExp
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
    PURE module function getLogSumExpMaxSeq_CK5(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG), intent(in)                :: maxArray
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogSumExpMaxSeq_CK4(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG), intent(in)                :: maxArray
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogSumExpMaxSeq_CK3(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG), intent(in)                :: maxArray
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogSumExpMaxSeq_CK2(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG), intent(in)                :: maxArray
        complex(CKG)                            :: logSumExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogSumExpMaxSeq_CK1(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: array(:)
        complex(CKG), intent(in)                :: maxArray
        complex(CKG)                            :: logSumExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogSumExpMaxSeq_RK5(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)   , intent(in)                :: maxArray
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogSumExpMaxSeq_RK4(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)   , intent(in)                :: maxArray
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogSumExpMaxSeq_RK3(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)   , intent(in)                :: maxArray
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogSumExpMaxSeq_RK2(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)   , intent(in)                :: maxArray
        real(RKG)                               :: logSumExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogSumExpMaxSeq_RK1(array, maxArray) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSeq_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: array(:)
        real(RKG)   , intent(in)                :: maxArray
        real(RKG)                               :: logSumExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLogSumExpMaxSel_CK5(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(selection_type), intent(in)                :: control
        complex(CKG)        , intent(in), contiguous    :: array(:)
        complex(CKG)        , intent(in)                :: maxArray
        complex(CKG)                                    :: logSumExp
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogSumExpMaxSel_CK4(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(selection_type), intent(in)                :: control
        complex(CKG)        , intent(in), contiguous    :: array(:)
        complex(CKG)        , intent(in)                :: maxArray
        complex(CKG)                                    :: logSumExp
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogSumExpMaxSel_CK3(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(selection_type), intent(in)                :: control
        complex(CKG)        , intent(in), contiguous    :: array(:)
        complex(CKG)        , intent(in)                :: maxArray
        complex(CKG)                                    :: logSumExp
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogSumExpMaxSel_CK2(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(selection_type), intent(in)                :: control
        complex(CKG)        , intent(in), contiguous    :: array(:)
        complex(CKG)        , intent(in)                :: maxArray
        complex(CKG)                                    :: logSumExp
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogSumExpMaxSel_CK1(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(selection_type), intent(in)                :: control
        complex(CKG)        , intent(in), contiguous    :: array(:)
        complex(CKG)        , intent(in)                :: maxArray
        complex(CKG)                                    :: logSumExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogSumExpMaxSel_RK5(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(selection_type), intent(in)                :: control
        real(RKG)           , intent(in), contiguous    :: array(:)
        real(RKG)           , intent(in)                :: maxArray
        real(RKG)                                       :: logSumExp
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogSumExpMaxSel_RK4(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(selection_type), intent(in)                :: control
        real(RKG)           , intent(in), contiguous    :: array(:)
        real(RKG)           , intent(in)                :: maxArray
        real(RKG)                                       :: logSumExp
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogSumExpMaxSel_RK3(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(selection_type), intent(in)                :: control
        real(RKG)           , intent(in), contiguous    :: array(:)
        real(RKG)           , intent(in)                :: maxArray
        real(RKG)                                       :: logSumExp
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogSumExpMaxSel_RK2(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(selection_type), intent(in)                :: control
        real(RKG)           , intent(in), contiguous    :: array(:)
        real(RKG)           , intent(in)                :: maxArray
        real(RKG)                                       :: logSumExp
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogSumExpMaxSel_RK1(array, maxArray, control) result(logSumExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExpMaxSel_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(selection_type), intent(in)                :: control
        real(RKG)           , intent(in), contiguous    :: array(:)
        real(RKG)           , intent(in)                :: maxArray
        real(RKG)                                       :: logSumExp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLogSumExp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathLogSumExp ! LCOV_EXCL_LINE
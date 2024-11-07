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
!>  This module contains the procedures and interfaces for computing the
!>  cumulative sum of the exponential of an array without undue numerical overflow.
!>
!>  \benchmarks
!>
!>  \benchmark{getCumPropExp_vs_setCumPropExp, The runtime performance of [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp) vs. [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)}
!>  \include{lineno} benchmark/pm_mathCumPropExp/getCumPropExp_vs_setCumPropExp/main.F90
!>  \compilefb{getCumPropExp_vs_setCumPropExp}
!>  \postprocb{getCumPropExp_vs_setCumPropExp}
!>  \include{lineno} benchmark/pm_mathCumPropExp/getCumPropExp_vs_setCumPropExp/main.py
!>  \visb{getCumPropExp_vs_setCumPropExp}
!>  \image html benchmark/pm_mathCumPropExp/getCumPropExp_vs_setCumPropExp/benchmark.getCumPropExp_vs_setCumPropExp.runtime.png width=1000
!>  \image html benchmark/pm_mathCumPropExp/getCumPropExp_vs_setCumPropExp/benchmark.getCumPropExp_vs_setCumPropExp.runtime.ratio.png width=1000
!>  \moralb{getCumPropExp_vs_setCumPropExp}
!>      -#  The procedures under the generic interface [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp) are functions while
!>          the procedures under the generic interface [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs significantly worse than the procedural interface.<br>
!>          However, the difference appears to diminish toward larger array sizes.<br>
!>
!>  \test
!>  [test_pm_mathCumPropExp](@ref test_pm_mathCumPropExp)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathCumPropExp

    use pm_control, only: control_type, sequence_type, sequence
    use pm_control, only: selection_type, selection
    use pm_array, only: direction_type, action_type
    use pm_array, only: backward, backward_type
    use pm_array, only: forward, forward_type
    use pm_array, only: reverse, reverse_type
    use pm_array, only: nothing, nothing_type
    use pm_kind, only: SK, IK, LK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathCumPropExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cumulative sum of the proportions of the exponential of the input array,
    !>  optionally in the backward direction and, optionally reverse the output cumulative sum upon return.
    !>
    !>  \details
    !>  The returned array is normalized such that all of its elements fall in the range \f$[0,1]\f$.<br>
    !>  All operations are performed while avoiding arithmetic overflow.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of type `real` of kind \RKALL
    !>                              whose cumulative proportional sum will have to be computed.
    !>  \param[in]  maxArray    :   The input scalar of the same type and kind as the input `array` representing
    !>                              the maximum value in `array` (i.e., `maxArray = maxval(array)`).
    !>  \param[in]  control     :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [sequence](@ref pm_control::sequence) or equivalently,
    !>                                          an object of type [sequence_type](@ref pm_control::sequence_type).<br>
    !>                                          Specifying this value forces the algorithm to skip runtime underflow checks.<br>
    !>                                          This means all exponentiation operations will be carried out for each element.<br>
    !>                                          Specifying this value can aid runtime efficiency when the divisions of none or very few
    !>                                          of the elements of `array` (for example, half or less) by `maxArray` causes underflow.<br>
    !>                                          In such cases, the potentially expensive runtime branching is avoided at the cost of performing a very few exponentiation operations.<br>
    !>                                          The typical cost of an if-branch is 7-20 CPU cycles on the contemporary architecture while exponentiation typically costs ~200 CPU cycles.<br>
    !>                                          See [the relevant benchmark here](#benchmark-setCumPropExp).<br>
    !>                                  <li>    the constant [selection](@ref pm_control::selection) or equivalently,
    !>                                          an object of type [selection_type](@ref pm_control::selection_type).<br>
    !>                                          Enabling this option can aid runtime efficiency when the division of a significant
    !>                                          number of elements of `array` (for example, half or more) by `maxArray` causes underflow.<br>
    !>                                          In such cases, the exponentiation is avoided if `control` = selection` leading to faster runtime
    !>                                          by avoiding exponentiation since it is highly expensive (on the order of ~200 CPU cycles).<br>
    !>                                          See [the relevant benchmark here](#benchmark-setCumPropExp).<br>
    !>                              </ol>
    !>                              (**optional**, default = [sequence](@ref pm_control::sequence))
    !>  \param[in]  direction   :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [forward](@ref pm_array::forward) or equivalently, an object of type [forward_type](@ref pm_array::forward_type),
    !>                                          implying that the output cumulative sum has be computed from the **first element** to the **last element** of the input `array`.<br>
    !>                                          even though the increments will still be written from the first element of `cumPropExp` to the last.<br>
    !>                                  <li>    the constant [backward](@ref pm_array::backward) or equivalently, an object of type [backward_type](@ref pm_array::backward_type),
    !>                                          implying that the output cumulative sum has be computed from the **last element** to the **first element** of the input `array`
    !>                                          even though the increments will still be written from the first element of `cumPropExp` to the last.<br>
    !>                              </ol>
    !>                              (**optional**, default = [sequence](@ref pm_control::sequence))
    !>  \param[in]  action      :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [nothing](@ref pm_array::nothing) or equivalently, an object of type [nothing_type](@ref pm_array::nothing_type),
    !>                                          implying no action to be performed on the elements of the output `cumPropExp` will have be reversed upon return.<br>
    !>                                  <li>    the constant [reverse](@ref pm_array::reverse) or equivalently, an object of type [reverse_type](@ref pm_array::reverse_type),
    !>                                          implying that the order of the elements of the output `cumPropExp` will have be reversed upon return,
    !>                                          such that its last element becomes the first.<br>
    !>                              </ol>
    !>                              (**optional**, default = [nothing](@ref pm_array::nothing))
    !>
    !>  \return
    !>  `cumPropExp`            :   The output array of the same size, shape, type, and kind as the input `array`
    !>                              containing the cumulative sum of proportions of `array` in the specified direction.
    !>
    !>  \interface{getCumPropExp}
    !>  \code{.F90}
    !>
    !>      use pm_mathCumPropExp, only: getCumPropExp, sequence, selection, forward, backward, nothing, reverse
    !>
    !>      cumPropExp(:) = getCumPropExp(array(:)          , maxArray = maxArray, direction = direction, action = action)
    !>      cumPropExp(:) = getCumPropExp(array(:), control , maxArray = maxArray, direction = direction, action = action)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `maxArray == maxval(array)` must hold for the corresponding arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  The functionalities of the procedures under this generic interface,
    !>  \code{.F90}
    !>      block
    !>          use pm_mathCumPropExp, only: getCumPropExp
    !>          cumPropExp = getCumPropExp(array)
    !>          cumPropExp = getCumPropExp(array, direction = backward)
    !>          cumPropExp = getCumPropExp(array, direction = backward, action = reversed)
    !>          cumPropExp = getCumPropExp(array, direction = forward, action = reversed)
    !>      end block
    !>      !
    !>  \endcode
    !>  are equivalent to the following lines respectively,
    !>  \code{.F90}
    !>      block
    !>          use pm_mathCumPropExp, only: getCumSum
    !>          cumPropExp = getCumSum(exp(array - maxval(array))); cumPropExp = cumPropExp / cumPropExp(size(array, maxArray, 1, IK))
    !>          cumPropExp = getCumSum(exp(array - maxval(array)), direction = backward); cumPropExp = cumPropExp / cumPropExp(size(array, maxArray, 1, IK))
    !>          cumPropExp = getCumSum(exp(array - maxval(array)), direction = backward, action = reversed); cumPropExp = cumPropExp / cumPropExp(1)
    !>          cumPropExp = getCumSum(exp(array - maxval(array)), direction = forward, action = reversed); cumPropExp = cumPropExp / cumPropExp(1)
    !>      end block
    !>      !
    !>  \endcode
    !>
    !>  \see
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [setCumSum](@ref pm_mathCumSum::setCumSum)<br>
    !>  [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp)<br>
    !>  [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)<br>
    !>
    !>  \example{getCumPropExp}
    !>  \include{lineno} example/pm_mathCumPropExp/getCumPropExp/main.F90
    !>  \compilef{getCumPropExp}
    !>  \output{getCumPropExp}
    !>  \include{lineno} example/pm_mathCumPropExp/getCumPropExp/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathCumPropExp](@ref test_pm_mathCumPropExp)
    !>
    !>  \final{getCumPropExp}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#if RK5_ENABLED
    PURE module function getCumPropExpDef_RK5(array, maxArray, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCumPropExpDef_RK4(array, maxArray, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCumPropExpDef_RK3(array, maxArray, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCumPropExpDef_RK2(array, maxArray, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCumPropExpDef_RK1(array, maxArray, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCumPropExpSel_RK5(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSel_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(selection_type)    , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCumPropExpSel_RK4(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSel_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(selection_type)    , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCumPropExpSel_RK3(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSel_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(selection_type)    , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCumPropExpSel_RK2(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSel_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(selection_type)    , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCumPropExpSel_RK1(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSel_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(selection_type)    , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCumPropExpSeq_RK5(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSeq_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(sequence_type)     , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCumPropExpSeq_RK4(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSeq_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(sequence_type)     , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCumPropExpSeq_RK3(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSeq_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(sequence_type)     , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCumPropExpSeq_RK2(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSeq_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(sequence_type)     , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCumPropExpSeq_RK1(array, maxArray, control, direction, action) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExpSeq_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: maxArray
        type(sequence_type)     , intent(in)                :: control
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumPropExp(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCumPropExp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the cumulative sum of the proportions of the exponential of the input array,
    !>  optionally in the backward direction and, optionally reverse the output cumulative sum upon return.
    !>
    !>  \details
    !>  The returned array is normalized such that all of its elements fall in the range \f$[0,1]\f$.<br>
    !>  All operations are performed while avoiding arithmetic overflow.
    !>
    !>  \param[out]     cumPropExp      :   The output array of the same size, shape, type, and kind as the input `array`
    !>                                      containing the cumulative sum of proportions of `array` in the specified direction.<br>
    !>                                      (**optional**, if missing, the result will be written to the input/output argument `array`.)
    !>  \param[inout]   array           :   The `contiguous` array of shape `(:)` of type `real` of kind \RKALL
    !>                                      whose cumulative proportional sum will have to be computed.<br>
    !>                                      <ul>
    !>                                          <li>    If `cumPropExp` is present, then `array` has `intent(in)`.<br>
    !>                                          <li>    If `cumPropExp` is missing, then `array` has `intent(inout)`.<br>
    !>                                                  On output, the contents of `array` will be completely overwritten by the computed `cumPropExp`.<br>
    !>                                      </ul>
    !>  \param[in]      maxArray        :   The input scalar of the same type and kind as the input `array` representing
    !>                                      the maximum value in `array` (i.e., `maxArray = maxval(array)`).
    !>  \param[in]      control         :   The input scalar object that can be,
    !>                                      <ol>
    !>                                          <li>    the constant [sequence](@ref pm_control::sequence) or equivalently,
    !>                                                  an object of type [sequence_type](@ref pm_control::sequence_type).<br>
    !>                                                  Specifying this value forces the algorithm to skip runtime underflow checks.<br>
    !>                                                  This means all exponentiation operations will be carried out for each element.<br>
    !>                                                  Specifying this value can aid runtime efficiency when the divisions of none or very few
    !>                                                  of the elements of `array` (for example, half or less) by `maxArray` causes underflow.<br>
    !>                                                  In such cases, the potentially expensive runtime branching is avoided at the cost of performing a very few exponentiation operations.<br>
    !>                                                  The typical cost of an if-branch is 7-20 CPU cycles on the contemporary architecture while exponentiation typically costs ~200 CPU cycles.<br>
    !>                                                  See [the relevant benchmark here](#benchmark-setCumPropExp).<br>
    !>                                          <li>    the constant [selection](@ref pm_control::selection) or equivalently,
    !>                                                  an object of type [selection_type](@ref pm_control::selection_type).<br>
    !>                                                  Enabling this option can aid runtime efficiency when the division of a significant
    !>                                                  number of elements of `array` (for example, half or more) by `maxArray` causes underflow.<br>
    !>                                                  In such cases, the exponentiation is avoided if `control` = selection` leading to faster runtime
    !>                                                  by avoiding exponentiation since it is highly expensive (on the order of ~200 CPU cycles).<br>
    !>                                                  See [the relevant benchmark here](#benchmark-setCumPropExp).<br>
    !>                                      </ol>
    !>  \param[in]      direction       :   The input scalar object that can be,
    !>                                      <ol>
    !>                                          <li>    the constant [forward](@ref pm_array::forward) or equivalently, an object of type [forward_type](@ref pm_array::forward_type),
    !>                                                  implying that the output cumulative sum has be computed from the **first element** to the **last element** of the input `array`.<br>
    !>                                                  even though the increments will still be written from the first element of `cumPropExp` to the last.<br>
    !>                                          <li>    the constant [backward](@ref pm_array::backward) or equivalently, an object of type [backward_type](@ref pm_array::backward_type),
    !>                                                  implying that the output cumulative sum has be computed from the **last element** to the **first element** of the input `array`
    !>                                                  even though the increments will still be written from the first element of `cumPropExp` to the last.<br>
    !>                                      </ol>
    !>                                      (**optional**, default = [sequence](@ref pm_control::sequence). It must be present **if and only if** the input argument `action` is also present.)
    !>  \param[in]      action          :   The input scalar object that can be,
    !>                                      <ol>
    !>                                          <li>    the constant [nothing](@ref pm_array::nothing) or equivalently, an object of type [nothing_type](@ref pm_array::nothing_type),
    !>                                                  implying no action to be performed on the elements of the output `cumPropExp` will have be reversed upon return.<br>
    !>                                          <li>    the constant [reverse](@ref pm_array::reverse) or equivalently, an object of type [reverse_type](@ref pm_array::reverse_type),
    !>                                                  implying that the order of the elements of the output `cumPropExp` will have be reversed upon return,
    !>                                                  such that its last element becomes the first.<br>
    !>                                      </ol>
    !>                                      (**optional**, default = [nothing](@ref pm_array::nothing). It must be present **if and only if** the input argument `direction` is also present.)
    !>
    !>  \interface{setCumPropExp}
    !>  \code{.F90}
    !>
    !>      use pm_mathCumPropExp, only: setCumPropExp
    !>      use pm_mathCumPropExp, only: selection, sequence, forward, backward, nothing, reverse
    !>
    !>      ! overwrite input array.
    !>
    !>      call setCumPropExp(array(:), maxArray, control)
    !>      call setCumPropExp(array(:), maxArray, control, direction, action)
    !>
    !>      ! write to new array.
    !>
    !>      call setCumPropExp(cumPropExp(:), array(:), maxArray, control)
    !>      call setCumPropExp(cumPropExp(:), array(:), maxArray, control, direction, action)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(array)` must hold for the corresponding arguments.<br>
    !>  The condition `maxArray == maxval(array)` must hold for the corresponding arguments.<br>
    !>  The condition `size(array) == size(cumPropExp)` must hold for the corresponding arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  The functionalities of the procedures under this generic interface,
    !>  \code{.F90}
    !>      block
    !>          use pm_mathCumPropExp, only: setCumPropExp
    !>          call setCumPropExp(cumPropExp, array, maxArray, control)
    !>          call setCumPropExp(cumPropExp, array, maxArray, control, direction = backward)
    !>          call setCumPropExp(cumPropExp, array, maxArray, control, direction = backward, action = reverse)
    !>          call setCumPropExp(cumPropExp, array, maxArray, control, direction = forward, action = reverse)
    !>      end block
    !>      !
    !>  \endcode
    !>  are equivalent to the following lines respectively,
    !>  \code{.F90}
    !>      block
    !>          use pm_mathCumPropExp, only: getCumSum
    !>          cumPropExp = getCumSum(exp(array - maxval(array))); cumPropExp = cumPropExp / cumPropExp(size(array, maxArray, 1, IK))
    !>          cumPropExp = getCumSum(exp(array - maxval(array)), direction = backward); cumPropExp = cumPropExp / cumPropExp(size(array, maxArray, 1, IK))
    !>          cumPropExp = getCumSum(exp(array - maxval(array)), direction = backward, action = reversed); cumPropExp = cumPropExp / cumPropExp(1)
    !>          cumPropExp = getCumSum(exp(array - maxval(array)), direction = forward, action = reversed); cumPropExp = cumPropExp / cumPropExp(1)
    !>      end block
    !>      !
    !>  \endcode
    !>
    !>  \see
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [setCumSum](@ref pm_mathCumSum::setCumSum)<br>
    !>  [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp)<br>
    !>  [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)<br>
    !>
    !>  \example{setCumPropExp}
    !>  \include{lineno} example/pm_mathCumPropExp/setCumPropExp/main.F90
    !>  \compilef{setCumPropExp}
    !>  \output{setCumPropExp}
    !>  \include{lineno} example/pm_mathCumPropExp/setCumPropExp/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setCumPropExp, The effects of `control` on runtime efficiency}
    !>      The following program compares the runtime performance of [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)
    !>      algorithm with and without checking for underflows.
    !>  \include{lineno} benchmark/pm_mathCumPropExp/setCumPropExp/main.F90
    !>  \compilefb{setCumPropExp}
    !>  \postprocb{setCumPropExp}
    !>  \include{lineno} benchmark/pm_mathCumPropExp/setCumPropExp/main.py
    !>  \visb{setCumPropExp}
    !>  \image html benchmark/pm_mathCumPropExp/setCumPropExp/benchmark.setCumPropExp.normal.png width=1000
    !>  \image html benchmark/pm_mathCumPropExp/setCumPropExp/benchmark.setCumPropExp.underflow.png width=1000
    !>  \moralb{setCumPropExp}
    !>      -#  If the input array has many (half the size of array or more) elements whose division by the `maxval(array)` causes underflow,
    !>          then setting `control = selection` when calling [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp) will likely result in a faster runtime.<br>
    !>          Conversely, if the divisions are not expected to cause any or too many underflows, then set `control = selection`
    !>          to improve cache coherence and runtime performance (at the expense of occasional expensive but redundant exponentiations).<br>
    !>      -#  If the input array size is less than 10-20 elements and [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp) is to be
    !>          called billions of times, then it would make sense to manually inline the procedure implementation in your code as
    !>          procedure call and processing of optional arguments will have a non-negligible performance overhead.
    !>
    !>  \test
    !>  [test_pm_mathCumPropExp](@ref test_pm_mathCumPropExp)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be expanded to include input arrays with `Weight`s.
    !>
    !>  \final{setCumPropExp}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

    ! seq, old, new, default

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqOldDefDef_RK5(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldDefDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqOldDefDef_RK4(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldDefDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqOldDefDef_RK3(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldDefDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqOldDefDef_RK2(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldDefDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqOldDefDef_RK1(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldDefDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqNewDefDef_RK5(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewDefDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqNewDefDef_RK4(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewDefDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqNewDefDef_RK3(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewDefDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqNewDefDef_RK2(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewDefDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqNewDefDef_RK1(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewDefDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! seq, old, forward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqOldForNon_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqOldForNon_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqOldForNon_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqOldForNon_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqOldForNon_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqOldForRev_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqOldForRev_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqOldForRev_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqOldForRev_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqOldForRev_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldForRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! seq, old, backward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacNon_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacNon_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacNon_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacNon_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacNon_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacRev_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacRev_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacRev_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacRev_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqOldBacRev_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqOldBacRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! seq, new, forward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqNewForNon_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqNewForNon_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqNewForNon_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqNewForNon_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqNewForNon_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqNewForRev_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqNewForRev_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqNewForRev_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqNewForRev_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqNewForRev_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewForRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! seq, new, backward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacNon_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacNon_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacNon_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacNon_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacNon_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacRev_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacRev_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacRev_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacRev_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSeqNewBacRev_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSeqNewBacRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(sequence_type) , intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! sel, old, new, default

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelOldDefDef_RK5(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldDefDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelOldDefDef_RK4(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldDefDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelOldDefDef_RK3(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldDefDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelOldDefDef_RK2(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldDefDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelOldDefDef_RK1(array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldDefDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelNewDefDef_RK5(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewDefDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelNewDefDef_RK4(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewDefDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelNewDefDef_RK3(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewDefDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelNewDefDef_RK2(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewDefDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelNewDefDef_RK1(cumPropExp, array, maxArray, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewDefDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! sel, old, forward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelOldForNon_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelOldForNon_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelOldForNon_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelOldForNon_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelOldForNon_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelOldForRev_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelOldForRev_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelOldForRev_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelOldForRev_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelOldForRev_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldForRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! sel, old, backward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelOldBacNon_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelOldBacNon_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelOldBacNon_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelOldBacNon_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelOldBacNon_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelOldBacRev_RK5(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelOldBacRev_RK4(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelOldBacRev_RK3(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelOldBacRev_RK2(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelOldBacRev_RK1(array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelOldBacRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! sel, new, forward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelNewForNon_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelNewForNon_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelNewForNon_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelNewForNon_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelNewForNon_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelNewForRev_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelNewForRev_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelNewForRev_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelNewForRev_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelNewForRev_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewForRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! sel, new, backward

    interface setCumPropExp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelNewBacNon_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelNewBacNon_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelNewBacNon_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelNewBacNon_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelNewBacNon_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCumPropExpSelNewBacRev_RK5(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCumPropExpSelNewBacRev_RK4(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCumPropExpSelNewBacRev_RK3(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCumPropExpSelNewBacRev_RK2(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCumPropExpSelNewBacRev_RK1(cumPropExp, array, maxArray, control, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumPropExpSelNewBacRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: cumPropExp(:)
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(in)                    :: maxArray
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
        type(selection_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Generate and return the normalized cumulative sum (i.e., Cumulative Density Function (CDF)) of the exponentials
    !>  of the input real vector robustly (without overflow or underflow). The last element of the returned vector is one.
    !>
    !>  \param[in]  array       :   The input `contiguous` vector of log-values whose log-sum-exp must be computed.
    !>  \param[in]  maxArray    :   The maximum of the input `array` argument (`maxArray = maxval(array)`).
    !>  \param[in]  lenArray    :   The length of the input array.
    !>
    !>  \return
    !>  `cumPropExp`            :   A real vector of the same length as the input array `array`.
    !>
    !>  \warning
    !>  This routine is only kept for backward compatibility and should not be used in production code.
    !>  Instead, use the procedures under the generic interface [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp).
    !>
    !>  \final
    !>
    !>  \test
    !>  [test_pm_mathCumPropExp](@ref test_pm_mathCumPropExp)
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    pure function getCumPropExp_RK(array, maxArray, lenArray) result(cumPropExp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumPropExp_RK
#endif
        use pm_kind, only: RKG => RK
        real(RKG)   , parameter     :: LOGTINY_RK = log(tiny(0._RKG))
        integer(IK) , intent(in)    :: lenArray
        real(RKG)   , intent(in)    :: array(lenArray)
        real(RKG)   , intent(in)    :: maxArray
        real(RKG)                   :: cumPropExp(lenArray)
        real(RKG)                   :: cumPropExpInv
        integer(IK)                 :: i
        cumPropExp(1) = array(1) - maxArray
        if (cumPropExp(1) < LOGTINY_RK) then
            cumPropExp(1) = 0._RKG
        else
            cumPropExp(1) = exp(cumPropExp(1))
        end if
        do i = 2, lenArray
            cumPropExp(i) = array(i) - maxArray
            if (cumPropExp(i) < LOGTINY_RK) then
                cumPropExp(i) = cumPropExp(i-1)
            else
                cumPropExp(i) = cumPropExp(i-1) + exp(cumPropExp(i))
            end if
        end do
        cumPropExpInv = 1._RKG / cumPropExp(lenArray)
        cumPropExp = cumPropExp * cumPropExpInv
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathCumPropExp ! LCOV_EXCL_LINE
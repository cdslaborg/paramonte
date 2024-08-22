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
!>  This module contains the procedures and interfaces for computing the cumulative sum of an array.
!>
!>  \benchmarks
!>
!>  \benchmark{getCumSum_vs_setCumSum, The runtime performance of [getCumSum](@ref pm_mathCumSum::getCumSum) vs. [setCumSum](@ref pm_mathCumSum::setCumSum)}
!>  \include{lineno} benchmark/pm_mathCumSum/getCumSum_vs_setCumSum/main.F90
!>  \compilefb{getCumSum_vs_setCumSum}
!>  \postprocb{getCumSum_vs_setCumSum}
!>  \include{lineno} benchmark/pm_mathCumSum/getCumSum_vs_setCumSum/main.py
!>  \visb{getCumSum_vs_setCumSum}
!>  \image html benchmark/pm_mathCumSum/getCumSum_vs_setCumSum/benchmark.getCumSum_vs_setCumSum.runtime.png width=1000
!>  \image html benchmark/pm_mathCumSum/getCumSum_vs_setCumSum/benchmark.getCumSum_vs_setCumSum.runtime.ratio.png width=1000
!>  \moralb{getCumSum_vs_setCumSum}
!>      -#  The procedures under the generic interface [getCumSum](@ref pm_mathCumSum::getCumSum) are functions while
!>          the procedures under the generic interface [setCumSum](@ref pm_mathCumSum::setCumSum) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs less efficiently than the subroutine interface.<br>
!>      -#  Furthermore, specifying the array bounds on the left-hand-side (LHS) assignment in the case of the functional interface (to avoid automatic reallocation)
!>          does not appear to enhance the performance of the functional interface in any meaningful way.<br>
!>          In other words, the compiler appears to be smart enough to not reallocate the LHS needlessly.<br>
!>
!>  \test
!>  [test_pm_mathCumSum](@ref test_pm_mathCumSum)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathCumSum

    use pm_array, only: direction_type, action_type
    use pm_array, only: backward, backward_type
    use pm_array, only: forward, forward_type
    use pm_array, only: reverse, reverse_type
    use pm_array, only: nothing, nothing_type
    use pm_kind, only: SK, IK, LK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathCumSum"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cumulative sum of the input array, optionally in the backward direction and,
    !>  optionally reverse the output cumulative sum array upon return.
    !>
    !>  \param[in]      array       :   The input `contiguous` array of shape `(:)` of either,
    !>                                  <ul>
    !>                                      <li>    type `integer` of kind \IKALL or, <br>
    !>                                      <li>    type `complex` of kind \CKALL or, <br>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ul>
    !>                                  whose cumulative sum will have to be computed.
    !>  \param[in]      direction   :   The input scalar object that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [forward](@ref pm_array::forward) or equivalently, an object of type [forward_type](@ref pm_array::forward_type),
    !>                                              implying that the output cumulative sum has be computed from the **first element** to the **last element** of the input `array`.<br>
    !>                                              even though the increments will still be written from the first element of `cumsum` to the last.<br>
    !>                                      <li>    the constant [backward](@ref pm_array::backward) or equivalently, an object of type [backward_type](@ref pm_array::backward_type),
    !>                                              implying that the output cumulative sum has be computed from the **last element** to the **first element** of the input `array`
    !>                                              even though the increments will still be written from the first element of `cumsum` to the last.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [forward](@ref pm_array::forward))
    !>  \param[in]      action      :   The input scalar object that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [nothing](@ref pm_array::nothing) or equivalently, an object of type [nothing_type](@ref pm_array::nothing_type),
    !>                                              implying no action to be performed on the elements of the output `cumsum`.<br>
    !>                                      <li>    the constant [reverse](@ref pm_array::reverse) or equivalently, an object of type [reverse_type](@ref pm_array::reverse_type),
    !>                                              implying that the order of the elements of the output `cumsum` will have be reversed upon return,
    !>                                              such that its last element becomes the first.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [nothing](@ref pm_array::nothing))
    !>
    !>  \return
    !>  `cumsum`                    :   The output array of the same size, shape, type, and kind as the input `array`
    !>                                  containing the cumulative sum of `array` in the specified direction.
    !>
    !>  \interface{getCumSum}
    !>  \code{.F90}
    !>
    !>      use pm_mathCumSum, only: getCumSum
    !>
    !>      cumsum(:) = getCumSum(array(:), direction = direction, action = action)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setCumSum](@ref pm_mathCumSum::setCumSum)<br>
    !>  [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp)<br>
    !>  [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)<br>
    !>
    !>  \example{getCumSum}
    !>  \include{lineno} example/pm_mathCumSum/getCumSum/main.F90
    !>  \compilef{getCumSum}
    !>  \output{getCumSum}
    !>  \include{lineno} example/pm_mathCumSum/getCumSum/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathCumSum](@ref test_pm_mathCumSum)
    !>
    !>  \todo
    !>  \pmed This generic interface can be expanded to include input arrays with `Weight`s.
    !>
    !>  \final{getCumSum}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface getCumSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCumSum_IK5(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        integer(IKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getCumSum_IK4(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        integer(IKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getCumSum_IK3(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        integer(IKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getCumSum_IK2(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        integer(IKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getCumSum_IK1(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        integer(IKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCumSum_CK5(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        complex(CKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCumSum_CK4(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        complex(CKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCumSum_CK3(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        complex(CKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCumSum_CK2(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        complex(CKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCumSum_CK1(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        complex(CKG)                                        :: cumsum(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCumSum_RK5(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumsum(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCumSum_RK4(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumsum(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCumSum_RK3(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumsum(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCumSum_RK2(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumsum(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCumSum_RK1(array, direction, action) result(cumsum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        class(direction_type)   , intent(in), optional      :: direction
        class(action_type)      , intent(in), optional      :: action
        real(RKG)                                           :: cumsum(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCumSum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the cumulative sum of the input `array`, optionally in the backward direction and optionally,
    !>  reverse the output cumulative sum array upon return.
    !>
    !>  \param[out]     cumsum      :   The output `contiguous` array of the same size, shape, type, and kind as the input `array`
    !>                                  containing the cumulative sum of `array` in the specified direction.<br>
    !>                                  (**optional**, if missing, the result will be written to the input/output argument `array`.)
    !>  \param[inout]   array       :   The `contiguous` array of shape `(:)` of either,
    !>                                  <ol>
    !>                                      <li>    type `integer` of kind \IKALL or, <br>
    !>                                      <li>    type `complex` of kind \CKALL or, <br>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  whose cumulative sum will have to be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `cumsum` is present, then `array` has `intent(in)`.<br>
    !>                                      <li>    If `cumsum` is missing, then `array` has `intent(inout)`.<br>
    !>                                              On output, the contents of `array` will be completely overwritten by the computed `cumsum`.<br>
    !>                                  </ol>
    !>  \param[in]      direction   :   The input scalar object that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [forward](@ref pm_array::forward) or equivalently, an object of type [forward_type](@ref pm_array::forward_type),
    !>                                              implying that the output cumulative sum has be computed from the **first element** to the **last element** of the input `array`.<br>
    !>                                              even though the increments will still be written from the first element of `cumsum` to the last.<br>
    !>                                      <li>    the constant [backward](@ref pm_array::backward) or equivalently, an object of type [backward_type](@ref pm_array::backward_type),
    !>                                              implying that the output cumulative sum has be computed from the **last element** to the **first element** of the input `array`
    !>                                              even though the increments will still be written from the first element of `cumsum` to the last.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `action` is also present.)
    !>  \param[in]      action      :   The input scalar object that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [nothing](@ref pm_array::nothing) or equivalently, an object of type [nothing_type](@ref pm_array::nothing_type),
    !>                                              implying no action to be performed on the elements of the output `cumsum`.<br>
    !>                                      <li>    the constant [reverse](@ref pm_array::reverse) or equivalently, an object of type [reverse_type](@ref pm_array::reverse_type),
    !>                                              implying that the order of the elements of the output `cumsum` will have be reversed upon return,
    !>                                              such that its last element becomes the first.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `direction` is also present.)
    !>
    !>  \interface{setCumSum}
    !>  \code{.F90}
    !>
    !>      use pm_mathCumSum, only: setCumSum
    !>
    !>      call setCumSum(array(:))
    !>      call setCumSum(cumsum(:), array(:))
    !>      call setCumSum(array(:), direction, action)
    !>      call setCumSum(cumsum(:), array(:), direction, action)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(array)` must hold for the corresponding arguments.<br>
    !>  The condition `size(array) == size(cumsum)` must hold for the corresponding arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \recursive
    !>
    !>  \remark
    !>  The procedures under this generic interface becomes `recursive` when `array` has `intent(inout)` (i.e., `cumsum` is missing).
    !>
    !>  \see
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [getCumPropExp](@ref pm_mathCumPropExp::getCumPropExp)<br>
    !>  [setCumPropExp](@ref pm_mathCumPropExp::setCumPropExp)<br>
    !>
    !>  \example{setCumSum}
    !>  \include{lineno} example/pm_mathCumSum/setCumSum/main.F90
    !>  \compilef{setCumSum}
    !>  \output{setCumSum}
    !>  \include{lineno} example/pm_mathCumSum/setCumSum/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathCumSum](@ref test_pm_mathCumSum)
    !>
    !>  \todo
    !>  \pmed This generic interface can be expanded to include input arrays with `Weight`s.
    !>
    !>  \final{setCumSum}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX

    ! old, new, default

    interface setCumSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_IK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_IK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_IK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_IK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_IK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_CK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_CK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_CK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_CK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_CK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_RK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_RK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_RK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_RK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumOldDefDef_RK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldDefDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_IK5(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_IK4(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_IK3(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_IK2(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_IK1(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_CK5(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_CK4(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_CK3(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_CK2(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_CK1(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_RK5(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_RK4(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_RK3(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_RK2(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumNewDefDef_RK1(cumsum, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewDefDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCumSum

    ! old, forward

    interface setCumSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_IK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_IK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_IK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_IK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_IK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_CK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_CK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_CK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_CK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_CK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_RK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_RK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_RK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_RK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumOldForNon_RK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_IK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_IK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_IK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_IK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_IK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_CK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_CK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_CK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_CK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_CK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_RK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_RK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_RK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_RK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumOldForRev_RK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldForRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCumSum

    ! old, backward

    interface setCumSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_IK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_IK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_IK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_IK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_IK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_CK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_CK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_CK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_CK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_CK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_RK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_RK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_RK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_RK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumOldBacNon_RK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_IK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_IK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_IK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_IK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_IK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_CK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_CK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_CK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_CK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_CK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_RK5(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_RK4(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_RK3(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_RK2(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumOldBacRev_RK1(array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumOldBacRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: array(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCumSum

    ! new, forward

    interface setCumSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_IK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_IK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_IK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_IK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_IK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_CK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_CK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_CK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_CK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_CK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_RK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_RK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_RK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_RK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumNewForNon_RK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_IK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_IK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_IK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_IK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_IK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_CK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_CK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_CK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_CK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_CK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_RK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_RK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_RK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_RK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumNewForRev_RK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewForRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(forward_type)  , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCumSum

    ! new, backward

    interface setCumSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_IK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_IK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_IK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_IK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_IK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_CK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_CK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_CK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_CK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_CK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_RK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_RK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_RK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_RK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumNewBacNon_RK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacNon_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(nothing_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_IK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_IK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_IK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_IK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_IK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)    , contiguous    :: array(:)
        integer(IKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_CK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_CK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_CK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_CK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_CK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: array(:)
        complex(CKG)        , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_RK5(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_RK4(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_RK3(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_RK2(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setCumSumNewBacRev_RK1(cumsum, array, direction, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCumSumNewBacRev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: array(:)
        real(RKG)           , intent(out)   , contiguous    :: cumsum(:)
        type(backward_type) , intent(in)                    :: direction
        type(reverse_type)  , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCumSum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathCumSum ! LCOV_EXCL_LINE
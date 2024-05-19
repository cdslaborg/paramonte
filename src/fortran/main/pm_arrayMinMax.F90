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
!>  This module contains procedures and generic interfaces for finding the minimum and maximum of two input **scalar** numbers through lexical comparison.
!>
!>  \details
!>  See [pm_mathMinMax](@ref pm_mathMinMax) for the equivalent operation on a pair of scalar values.<br>
!>
!>  \see
!>  [pm_swap](@ref pm_swap)<br>
!>  [pm_mathDivMul](@ref pm_mathDivMul)<br>
!>  [pm_mathSubAdd](@ref pm_mathSubAdd)<br>
!>  [pm_mathMinMax](@ref pm_mathMinMax) (for lexical minimum/maximum of a pair of **scalar** values)<br>
!>
!>  \final
!>
!>  \test
!>  [test_pm_arrayMinMax](@ref test_pm_arrayMinMax)
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayMinMax

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayMinMax"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an array of size two containing the minimum and maximum values of a sequence of values in the first and second elements, respectively.<br>
    !>
    !>  \details
    !>  Note that complex values are compared lexicographically.<br>
    !>  However, when the input sequence is a scalar string, the individual characters comprise the elements/values that are compared against each other.<br>
    !>  See [pm_complexCompareLex](@ref pm_complexCompareLex) and [pm_container](@ref pm_container) for more information on the relevant lexical comparison operators.<br>
    !>  Also, note that a logical `.false.` is less than `.true.` in this module.<br>
    !>
    !>  \param[in]  array   :   The input scalar sequence of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter,<br>
    !>                          </ol>
    !>                          or an array of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter,<br>
    !>                              <li>    type `integer` of kind \IKALL,<br>
    !>                              <li>    type `logical` of kind \LKALL,<br>
    !>                              <li>    type `complex` of kind \CKALL,<br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt) of kind \SKALL or<br>
    !>                              <li>    type [css_type](@ref pm_container::css_type) of default kind \SK,<br>
    !>                          </ol>
    !>                          containing the values whose lexicographic minimum and maximum values are to be found.<br>
    !>
    !>  \return
    !>  `minMaxVal`         :   The output scalar of,
    !>                          <ol>
    !>                              <li>    type `character` of length `2` of the same kind as the input **scalar** `array` rank `0` of type `character`,<br>
    !>                          </ol>
    !>                          or otherwise, vector of shape `(1:2)` of the same type and kind (and length-type parameter) as the input `array` of rank `1`,
    !>                          containing `minval(array)` in the first element and `maxval(array) in the second element.<br>
    !>                          <ol>
    !>                              <li>    If the input `array` is a scalar of type `character` of arbitrary length-type parameter, then `minMaxVal(1)`
    !>                                      will be set to the largest value for the corresponding `character` kind in the processor collating sequence.<br>
    !>                              <li>    If the input `array` is a vector of size zero, then `minMaxVal(1)` will be set to the largest value for the kind of `minMaxVal(1)`.<br>
    !>                              <li>    If the input `array` is a vector of containers, then the `val` component of `minMaxVal(1)` will remain **unallocated** on output.<br>
    !>                          </ol>
    !>                          <ol>
    !>                              <li>    If the input `array` is a scalar of type `character` of arbitrary length-type parameter, then `minMaxVal(2)`
    !>                                      will be set to the smallest value for the corresponding `character` kind in the processor collating sequence.<br>
    !>                              <li>    If the input `array` is a vector of size zero, then `minMaxVal(2)` will be set to the smallest value for the kind of `minMaxVal(2)`.<br>
    !>                              <li>    If the input `array` is a vector of containers, then the `val` component of `minMaxVal(2)` will remain **unallocated** on output.<br>
    !>                          </ol>
    !>
    !>  \interface{getMinMaxVal}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMinMax, only: getMinMaxVal
    !>
    !>      minMaxVal(1:2) = getMinMaxVal(array) ! on output, minMaxVal(1) <= minMaxVal(2) holds.
    !>      !
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \note
    !>  The procedures under this generic interface are particularly useful in combination
    !>  with procedures of [pm_arrayMembership](@ref pm_arrayMembership) and similar interfaces.<br>
    !>
    !>  \see
    !>  [setMinMaxVal](@ref pm_arrayMinMax::setMinMaxVal)<br>
    !>  [pm_arrayMembership](@ref pm_arrayMembership)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{getMinMaxVal}
    !>  \include{lineno} example/pm_arrayMinMax/getMinMaxVal/main.F90
    !>  \compilef{getMinMaxVal}
    !>  \output{getMinMaxVal}
    !>  \include{lineno} example/pm_arrayMinMax/getMinMaxVal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMinMax](@ref test_pm_arrayMinMax)
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to include the possibility of passing user-defined custom comparison.
    !>
    !>  \final{getMinMaxVal}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface getMinMaxVal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getMinMaxVal_D0_SK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(2,SKG)                                        :: minMaxVal
    end function
#endif

#if SK4_ENABLED
    pure module function getMinMaxVal_D0_SK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(2,SKG)                                        :: minMaxVal
    end function
#endif

#if SK3_ENABLED
    pure module function getMinMaxVal_D0_SK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(2,SKG)                                        :: minMaxVal
    end function
#endif

#if SK2_ENABLED
    pure module function getMinMaxVal_D0_SK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(2,SKG)                                        :: minMaxVal
    end function
#endif

#if SK1_ENABLED
    pure module function getMinMaxVal_D0_SK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(2,SKG)                                        :: minMaxVal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getMinMaxVal_D1_SK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: minMaxVal(2)
    end function
#endif

#if SK4_ENABLED
    pure module function getMinMaxVal_D1_SK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: minMaxVal(2)
    end function
#endif

#if SK3_ENABLED
    pure module function getMinMaxVal_D1_SK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: minMaxVal(2)
    end function
#endif

#if SK2_ENABLED
    pure module function getMinMaxVal_D1_SK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: minMaxVal(2)
    end function
#endif

#if SK1_ENABLED
    pure module function getMinMaxVal_D1_SK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: minMaxVal(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getMinMaxVal_D1_IK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: minMaxVal(2)
    end function
#endif

#if IK4_ENABLED
    pure module function getMinMaxVal_D1_IK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: minMaxVal(2)
    end function
#endif

#if IK3_ENABLED
    pure module function getMinMaxVal_D1_IK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: minMaxVal(2)
    end function
#endif

#if IK2_ENABLED
    pure module function getMinMaxVal_D1_IK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: minMaxVal(2)
    end function
#endif

#if IK1_ENABLED
    pure module function getMinMaxVal_D1_IK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: minMaxVal(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function getMinMaxVal_D1_LK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: minMaxVal(2)
    end function
#endif

#if LK4_ENABLED
    pure module function getMinMaxVal_D1_LK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: minMaxVal(2)
    end function
#endif

#if LK3_ENABLED
    pure module function getMinMaxVal_D1_LK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: minMaxVal(2)
    end function
#endif

#if LK2_ENABLED
    pure module function getMinMaxVal_D1_LK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: minMaxVal(2)
    end function
#endif

#if LK1_ENABLED
    pure module function getMinMaxVal_D1_LK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: minMaxVal(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getMinMaxVal_D1_CK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: minMaxVal(2)
    end function
#endif

#if CK4_ENABLED
    pure module function getMinMaxVal_D1_CK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: minMaxVal(2)
    end function
#endif

#if CK3_ENABLED
    pure module function getMinMaxVal_D1_CK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: minMaxVal(2)
    end function
#endif

#if CK2_ENABLED
    pure module function getMinMaxVal_D1_CK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: minMaxVal(2)
    end function
#endif

#if CK1_ENABLED
    pure module function getMinMaxVal_D1_CK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: minMaxVal(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getMinMaxVal_D1_RK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: minMaxVal(2)
    end function
#endif

#if RK4_ENABLED
    pure module function getMinMaxVal_D1_RK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: minMaxVal(2)
    end function
#endif

#if RK3_ENABLED
    pure module function getMinMaxVal_D1_RK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: minMaxVal(2)
    end function
#endif

#if RK2_ENABLED
    pure module function getMinMaxVal_D1_RK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: minMaxVal(2)
    end function
#endif

#if RK1_ENABLED
    pure module function getMinMaxVal_D1_RK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: minMaxVal(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getMinMaxVal_D1_PSSK5(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: minMaxVal(2)
    end function
#endif

#if SK4_ENABLED
    module function getMinMaxVal_D1_PSSK4(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: minMaxVal(2)
    end function
#endif

#if SK3_ENABLED
    module function getMinMaxVal_D1_PSSK3(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: minMaxVal(2)
    end function
#endif

#if SK2_ENABLED
    module function getMinMaxVal_D1_PSSK2(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: minMaxVal(2)
    end function
#endif

#if SK1_ENABLED
    module function getMinMaxVal_D1_PSSK1(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: minMaxVal(2)
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getMinMaxVal_D1_BSSK(array) result(minMaxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxVal_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        type(css_type)                                          :: minMaxVal(2)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getMinMaxVal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the minimum and maximum values of the input sequence.<br>
    !>
    !>  \brief
    !>  Note that complex values are compared lexicographically.<br>
    !>  However, when the input sequence is a scalar string, the individual characters comprise the elements/values that are compared against each other.<br>
    !>  See [pm_complexCompareLex](@ref pm_complexCompareLex) and [pm_container](@ref pm_container) for more information on the relevant lexical comparison operators.<br>
    !>  Also, note that a logical `.false.` is less than `.true.` in this module.<br>
    !>
    !>  \param[in]  array   :   The input scalar sequence of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter,<br>
    !>                          </ol>
    !>                          or an array of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter,<br>
    !>                              <li>    type `integer` of kind \IKALL,<br>
    !>                              <li>    type `logical` of kind \LKALL,<br>
    !>                              <li>    type `complex` of kind \CKALL,<br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt) of kind \SKALL or<br>
    !>                              <li>    type [css_type](@ref pm_container::css_type) of default kind \SK,<br>
    !>                          </ol>
    !>                          containing the values whose lexicographic minimum and maximum values are to be found.<br>
    !>  \param[out] vmin    :   The output scalar of,
    !>                          <ol>
    !>                              <li>    type `character` of length `1` of the same kind as the input **scalar** `array` rank `0` of type `character`,<br>
    !>                          </ol>
    !>                          or otherwise, of the same type and kind (and length-type parameter) as the input `array` of rank `1`,
    !>                          containing the minimum value in the sequence.<br>
    !>                          <ol>
    !>                              <li>    If the input `array` is a scalar of type `character` of arbitrary length-type parameter, then `vmin` will be set to the largest value for the corresponding `character` kind in the processor collating sequence.<br>
    !>                              <li>    If the input `array` is a vector of size zero, then `vmin` will be set to the largest value for the kind of `vmin`.<br>
    !>                              <li>    If the input `array` is a vector of containers, then the `val` component of `vmin` will remain **unallocated** on output.<br>
    !>                          </ol>
    !>  \param[out] vmax    :   The output scalar of,
    !>                          <ol>
    !>                              <li>    type `character` of length `1` of the same kind as the input **scalar** `array` rank `0` of type `character`,<br>
    !>                          </ol>
    !>                          or otherwise, of the same type and kind (and length-type parameter) as the input `array` of rank `1`,
    !>                          containing the maximum value in the sequence.<br>
    !>                          <ol>
    !>                              <li>    If the input `array` is a scalar of type `character` of arbitrary length-type parameter, then `vmax` will be set to the smallest value for the corresponding `character` kind in the processor collating sequence.<br>
    !>                              <li>    If the input `array` is a vector of size zero, then `vmax` will be set to the smallest value for the kind of `vmax`.<br>
    !>                              <li>    If the input `array` is a vector of containers, then the `val` component of `vmax` will remain **unallocated** on output.<br>
    !>                          </ol>
    !>
    !>  \interface{setMinMaxVal}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMinMax, only: setMinMaxVal
    !>
    !>      call setMinMaxVal(array, vmin, vmax)
    !>      !
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \note
    !>  The procedures under this generic interface are particularly useful in combination
    !>  with procedures of [pm_arrayMembership](@ref pm_arrayMembership) and similar interfaces.<br>
    !>
    !>  \see
    !>  [setMinMaxVal](@ref pm_arrayMinMax::setMinMaxVal)<br>
    !>  [pm_arrayMembership](@ref pm_arrayMembership)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{setMinMaxVal}
    !>  \include{lineno} example/pm_arrayMinMax/setMinMaxVal/main.F90
    !>  \compilef{setMinMaxVal}
    !>  \output{setMinMaxVal}
    !>  \include{lineno} example/pm_arrayMinMax/setMinMaxVal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMinMax](@ref test_pm_arrayMinMax)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.8.0 20221119}, \ifx{2023.0.0 20221201}
    !>  \desc
    !>  The \ifort and \ifx as of the specified versions above cannot find the minimum and maximum of a constant empty array of type `character`.<br>
    !>  This has prevented the definition of the constants `COL_SEQ_BEG` and `COL_SEQ_END` in the module [pm_str](@ref pm_str) which are also needed in this module.<br>
    !>  For example, the following code yields an Internal Compiler Error:
    !>  \code{.F90}
    !>      character(1), parameter :: ARRAY_EMPTY(0) = [character(1)::]
    !>      character(*), parameter :: COL_SEQ_BEG = minval(ARRAY_EMPTY,1)
    !>      end
    !>  \endcode
    !>  \verbatim
    !>      catastrophic error: **Internal compiler error: segmentation violation signal raised**
    !>      Please report this error along with the circumstances in which it occurred in a Software Problem Report.
    !>      Note: File and line given may not be explicit cause of this error.
    !>  \endverbatim
    !>  \remedy
    !>  For now, the empty arrays are constructed dynamically within the routines of this module.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The dynamic `character` allocation due to the bug described above must be eliminated
    !>  within the implementation of this module as soon as the Intel compilers bugs are resolved.
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to include the possibility of passing user-defined custom comparison.
    !>
    !>  \final{setMinMaxVal}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

    interface setMinMaxVal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setMinMaxVal_D0_SK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                        , intent(in)                    :: array
        character(1,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setMinMaxVal_D0_SK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                        , intent(in)                    :: array
        character(1,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setMinMaxVal_D0_SK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                        , intent(in)                    :: array
        character(1,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setMinMaxVal_D0_SK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                        , intent(in)                    :: array
        character(1,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setMinMaxVal_D0_SK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                        , intent(in)                    :: array
        character(1,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setMinMaxVal_D1_SK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setMinMaxVal_D1_SK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setMinMaxVal_D1_SK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setMinMaxVal_D1_SK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setMinMaxVal_D1_SK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)                        , intent(out)                   :: vmin, vmax
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setMinMaxVal_D1_IK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setMinMaxVal_D1_IK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setMinMaxVal_D1_IK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setMinMaxVal_D1_IK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setMinMaxVal_D1_IK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setMinMaxVal_D1_LK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setMinMaxVal_D1_LK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setMinMaxVal_D1_LK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setMinMaxVal_D1_LK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setMinMaxVal_D1_LK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setMinMaxVal_D1_CK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setMinMaxVal_D1_CK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setMinMaxVal_D1_CK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setMinMaxVal_D1_CK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setMinMaxVal_D1_CK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                            , intent(out)                   :: vmin, vmax
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setMinMaxVal_D1_RK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                               , intent(in)    , contiguous    :: array(:)
        real(RKG)                               , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setMinMaxVal_D1_RK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                               , intent(in)    , contiguous    :: array(:)
        real(RKG)                               , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setMinMaxVal_D1_RK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                               , intent(in)    , contiguous    :: array(:)
        real(RKG)                               , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setMinMaxVal_D1_RK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                               , intent(in)    , contiguous    :: array(:)
        real(RKG)                               , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setMinMaxVal_D1_RK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                               , intent(in)    , contiguous    :: array(:)
        real(RKG)                               , intent(out)                   :: vmin, vmax
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setMinMaxVal_D1_PSSK5(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))                      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                      , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setMinMaxVal_D1_PSSK4(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))                      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                      , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setMinMaxVal_D1_PSSK3(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))                      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                      , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setMinMaxVal_D1_PSSK2(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))                      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                      , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setMinMaxVal_D1_PSSK1(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))                      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                      , intent(out)                   :: vmin, vmax
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setMinMaxVal_D1_BSSK(array, vmin, vmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxVal_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)                          , intent(in)    , contiguous    :: array(:)
        type(css_type)                          , intent(out)                   :: vmin, vmax
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setMinMaxVal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayMinMax ! LCOV_EXCL_LINE
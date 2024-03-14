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
!>  This module contains procedures and generic interfaces for finding the minimum and maximum of two input **scalar** values through lexical comparison.
!>
!>  \details
!>  See [pm_arrayMinMax](@ref pm_arrayMinMax) for the equivalent operation on a **sequence** of values.<br>
!>
!>  \see
!>  [pm_swap](@ref pm_swap)<br>
!>  [pm_mathDivMul](@ref pm_mathDivMul)<br>
!>  [pm_mathSubAdd](@ref pm_mathSubAdd)<br>
!>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
!>
!>  \finmain
!>
!>  \test
!>  [test_pm_mathMinMax](@ref test_pm_mathMinMax)
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathMinMax

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathMinMax"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an array of size two containing the minimum and maximum of the two input values in the first and second elements, respectively.<br>
    !>
    !>  \brief
    !>  Note that strings and complex values are compared lexicographically.<br>
    !>  See [pm_complexCompareLex](@ref pm_complexCompareLex) and [pm_container](@ref pm_container) for more information on the relevant lexical comparison operators.<br>
    !>  Also, note that a logical `.false.` is less than `.true.` in this module.<br>
    !>
    !>  \param[in]  a       :   The input scalar of the same type and kind as the output `minMax`.<br>
    !>                          (**optional**, it must be present <b>if and only if</b> the input argument `b` is also present.)
    !>  \param[in]  b       :   The input scalar of the same type and kind as the output `minMax`.<br>
    !>                          (**optional**, it must be present <b>if and only if</b> the input argument `a` is also present.)
    !>  \param[in]  pair    :   The input vector of length `2` of the same type and kind as the output `minMax`.<br>
    !>                          (**optional**, it must be present <b>if and only if</b> the input arguments `a` and `b` are missing.)
    !>
    !>  \return
    !>  `minMax`            :   The output vector of shape `(2)` of either<br>
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter, or<br>
    !>                              <li>    type `integer` of kind \IKALL or <br>
    !>                              <li>    type `logical` of kind \LKALL or <br>
    !>                              <li>    type `complex` of kind \CKALL or <br>
    !>                              <li>    type `real` of kind \RKALL.<br>
    !>                          </ol>
    !>                          containing `min(a,b)` or `minval(pair)` in the first element and `max(a,b)` or `maxval(pair)` in the second element.<br>
    !>
    !>  \interface{getMinMax}
    !>  \code{.F90}
    !>
    !>      use pm_mathMinMax, only: getMinMax
    !>
    !>      minMax(1:2) = getMinMax(a, b) ! on output, minMax(1) <= minMax(2) holds.
    !>      minMax(1:2) = getMinMax(pair(1:2)) ! on output, minMax(1) <= minMax(2) holds.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  If the two input arguments `a` and `b` are of type `character`, then they are assumed to have the same length.<br>
    !>  Otherwise, the output `minMax` will have a length type parameter that is the maximum of the two input values.<br>
    !>  This means that array padding with blanks will occur for the value of the shorter string.<br>
    !>  For clarifications, see the examples below.<br>
    !>
    !>  \pure
    !>
    !>  \note
    !>  The procedures under this generic interface are particularly useful in combination
    !>  with [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp) and similar interfaces.<br>
    !>
    !>  \see
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>
    !>  \example{getMinMax}
    !>  \include{lineno} example/pm_mathMinMax/getMinMax/main.F90
    !>  \compilef{getMinMax}
    !>  \output{getMinMax}
    !>  \include{lineno} example/pm_mathMinMax/getMinMax/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathMinMax](@ref test_pm_mathMinMax)
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to include the possibility of passing user-defined custom comparison.
    !>
    !>  \finmain{getMinMax}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface getMinMax

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getMinMaxIndi_SK5(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(in)                :: a, b
        character(max(len(a,IK),len(b,IK)),SKC)                             :: minMax(2)
    end function
#endif

#if SK4_ENABLED
    pure module function getMinMaxIndi_SK4(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(in)                :: a, b
        character(max(len(a,IK),len(b,IK)),SKC)                             :: minMax(2)
    end function
#endif

#if SK3_ENABLED
    pure module function getMinMaxIndi_SK3(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(in)                :: a, b
        character(max(len(a,IK),len(b,IK)),SKC)                             :: minMax(2)
    end function
#endif

#if SK2_ENABLED
    pure module function getMinMaxIndi_SK2(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(in)                :: a, b
        character(max(len(a,IK),len(b,IK)),SKC)                             :: minMax(2)
    end function
#endif

#if SK1_ENABLED
    pure module function getMinMaxIndi_SK1(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(in)                :: a, b
        character(max(len(a,IK),len(b,IK)),SKC)                             :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getMinMaxIndi_IK5(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(in)                :: a, b
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK4_ENABLED
    pure module function getMinMaxIndi_IK4(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(in)                :: a, b
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK3_ENABLED
    pure module function getMinMaxIndi_IK3(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(in)                :: a, b
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK2_ENABLED
    pure module function getMinMaxIndi_IK2(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(in)                :: a, b
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK1_ENABLED
    pure module function getMinMaxIndi_IK1(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(in)                :: a, b
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function getMinMaxIndi_LK5(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(in)                :: a, b
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK4_ENABLED
    pure module function getMinMaxIndi_LK4(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(in)                :: a, b
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK3_ENABLED
    pure module function getMinMaxIndi_LK3(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(in)                :: a, b
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK2_ENABLED
    pure module function getMinMaxIndi_LK2(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(in)                :: a, b
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK1_ENABLED
    pure module function getMinMaxIndi_LK1(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(in)                :: a, b
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getMinMaxIndi_CK5(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(in)                :: a, b
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK4_ENABLED
    pure module function getMinMaxIndi_CK4(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(in)                :: a, b
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK3_ENABLED
    pure module function getMinMaxIndi_CK3(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(in)                :: a, b
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK2_ENABLED
    pure module function getMinMaxIndi_CK2(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(in)                :: a, b
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK1_ENABLED
    pure module function getMinMaxIndi_CK1(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(in)                :: a, b
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getMinMaxIndi_RK5(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(in)                :: a, b
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK4_ENABLED
    pure module function getMinMaxIndi_RK4(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(in)                :: a, b
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK3_ENABLED
    pure module function getMinMaxIndi_RK3(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(in)                :: a, b
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK2_ENABLED
    pure module function getMinMaxIndi_RK2(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(in)                :: a, b
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK1_ENABLED
    pure module function getMinMaxIndi_RK1(a, b) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxIndi_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(in)                :: a, b
        real(RKC)                                                           :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getMinMaxPair_SK5(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(in)                :: pair(2)
        character(len(pair,IK),SKC)                                         :: minMax(2)
    end function
#endif

#if SK4_ENABLED
    pure module function getMinMaxPair_SK4(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(in)                :: pair(2)
        character(len(pair,IK),SKC)                                         :: minMax(2)
    end function
#endif

#if SK3_ENABLED
    pure module function getMinMaxPair_SK3(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(in)                :: pair(2)
        character(len(pair,IK),SKC)                                         :: minMax(2)
    end function
#endif

#if SK2_ENABLED
    pure module function getMinMaxPair_SK2(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(in)                :: pair(2)
        character(len(pair,IK),SKC)                                         :: minMax(2)
    end function
#endif

#if SK1_ENABLED
    pure module function getMinMaxPair_SK1(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(in)                :: pair(2)
        character(len(pair,IK),SKC)                                         :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getMinMaxPair_IK5(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(in)                :: pair(2)
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK4_ENABLED
    pure module function getMinMaxPair_IK4(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(in)                :: pair(2)
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK3_ENABLED
    pure module function getMinMaxPair_IK3(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(in)                :: pair(2)
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK2_ENABLED
    pure module function getMinMaxPair_IK2(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(in)                :: pair(2)
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

#if IK1_ENABLED
    pure module function getMinMaxPair_IK1(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(in)                :: pair(2)
        integer(IKC)                                                        :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function getMinMaxPair_LK5(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(in)                :: pair(2)
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK4_ENABLED
    pure module function getMinMaxPair_LK4(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(in)                :: pair(2)
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK3_ENABLED
    pure module function getMinMaxPair_LK3(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(in)                :: pair(2)
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK2_ENABLED
    pure module function getMinMaxPair_LK2(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(in)                :: pair(2)
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

#if LK1_ENABLED
    pure module function getMinMaxPair_LK1(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(in)                :: pair(2)
        logical(LKC)                                                        :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getMinMaxPair_CK5(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(in)                :: pair(2)
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK4_ENABLED
    pure module function getMinMaxPair_CK4(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(in)                :: pair(2)
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK3_ENABLED
    pure module function getMinMaxPair_CK3(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(in)                :: pair(2)
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK2_ENABLED
    pure module function getMinMaxPair_CK2(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(in)                :: pair(2)
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

#if CK1_ENABLED
    pure module function getMinMaxPair_CK1(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(in)                :: pair(2)
        complex(CKC)                                                        :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getMinMaxPair_RK5(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(in)                :: pair(2)
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK4_ENABLED
    pure module function getMinMaxPair_RK4(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(in)                :: pair(2)
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK3_ENABLED
    pure module function getMinMaxPair_RK3(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(in)                :: pair(2)
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK2_ENABLED
    pure module function getMinMaxPair_RK2(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(in)                :: pair(2)
        real(RKC)                                                           :: minMax(2)
    end function
#endif

#if RK1_ENABLED
    pure module function getMinMaxPair_RK1(pair) result(minMax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinMaxPair_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(in)                :: pair(2)
        real(RKC)                                                           :: minMax(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getMinMax

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the minimum and maximum of the two input scalar values `a` and `b` in `a` and `b`, respectively, on output.<br>
    !>
    !>  \brief
    !>  Note that strings and complex values are compared lexicographically. See [pm_complexCompareLex](@ref pm_complexCompareLex)
    !>  and [pm_container](@ref pm_container) for more information on the relevant lexical comparison operators.<br>
    !>  Also, note that a logical `.false.` is less than `.true.` in this module.<br>
    !>
    !>  \param[inout]   a       :   The input/output scalar, or array of the same rank, shape, and size other array-like arguments, of either <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter, or<br>
    !>                                  <li>    type `integer` of kind \IKALL or <br>
    !>                                  <li>    type `logical` of kind \LKALL or <br>
    !>                                  <li>    type `complex` of kind \CKALL or <br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ol>
    !>                              On output, its value will be overwritten with `min(a,b)`.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> the argument `b` is also present.)
    !>  \param[inout]   b       :   The input/output scalar, or array of the same rank, shape, and size other array-like arguments, of the same type and kind as the input `a`.<br>
    !>                              On output, its value will be overwritten with `max(a,b)`.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> the argument `a` is also present.)
    !>  \param[inout]   pair    :   The input/output vector of size `2` of either, <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter, or<br>
    !>                                  <li>    type `integer` of kind \IKALL or <br>
    !>                                  <li>    type `logical` of kind \LKALL or <br>
    !>                                  <li>    type `complex` of kind \CKALL or <br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ol>
    !>                              On output, its value will be overwritten with `[min(a,b), max(a,b)]`.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> the arguments `a` and `b` are missing.)
    !>
    !>  \interface{setMinMax}
    !>  \code{.F90}
    !>
    !>      use pm_mathMinMax, only: setMinMax
    !>
    !>      call setMinMax(a, b) ! on output, a <= b holds.
    !>      call setMinMax(pair) ! on output, pair(1) <= pair(2) holds.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  If the two input arguments `a` and `b` are of type `character`, then they are assumed to have the same length.<br>
    !>  Otherwise, the output `minMax` will have a length type parameter that is the maximum of the two input values.<br>
    !>  This means that array padding with blanks will occur for the value of the shorter string.<br>
    !>  For clarifications, see the examples below.<br>
    !>
    !>  \warning
    !>  The rank and size of the two input arguments must be the same.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are not `elemental` when the input argument `pair` is present.<br>
    !>
    !>  \note
    !>  The procedures under this generic interface are particularly useful in combination
    !>  with [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp) and similar interfaces.<br>
    !>
    !>  \see
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>
    !>  \example{setMinMax}
    !>  \include{lineno} example/pm_mathMinMax/setMinMax/main.F90
    !>  \compilef{setMinMax}
    !>  \output{setMinMax}
    !>  \include{lineno} example/pm_mathMinMax/setMinMax/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathMinMax](@ref test_pm_mathMinMax)
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to include the possibility of passing user-defined custom comparison.
    !>
    !>  \finmain{setMinMax}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface setMinMax

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module subroutine setMinMaxIndi_SK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(inout)             :: a, b
    end subroutine
#endif

#if SK4_ENABLED
    pure elemental module subroutine setMinMaxIndi_SK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(inout)             :: a, b
    end subroutine
#endif

#if SK3_ENABLED
    pure elemental module subroutine setMinMaxIndi_SK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(inout)             :: a, b
    end subroutine
#endif

#if SK2_ENABLED
    pure elemental module subroutine setMinMaxIndi_SK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(inout)             :: a, b
    end subroutine
#endif

#if SK1_ENABLED
    pure elemental module subroutine setMinMaxIndi_SK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(inout)             :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module subroutine setMinMaxIndi_IK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if IK4_ENABLED
    pure elemental module subroutine setMinMaxIndi_IK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if IK3_ENABLED
    pure elemental module subroutine setMinMaxIndi_IK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if IK2_ENABLED
    pure elemental module subroutine setMinMaxIndi_IK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if IK1_ENABLED
    pure elemental module subroutine setMinMaxIndi_IK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module subroutine setMinMaxIndi_LK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if LK4_ENABLED
    pure elemental module subroutine setMinMaxIndi_LK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if LK3_ENABLED
    pure elemental module subroutine setMinMaxIndi_LK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if LK2_ENABLED
    pure elemental module subroutine setMinMaxIndi_LK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if LK1_ENABLED
    pure elemental module subroutine setMinMaxIndi_LK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module subroutine setMinMaxIndi_CK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if CK4_ENABLED
    pure elemental module subroutine setMinMaxIndi_CK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if CK3_ENABLED
    pure elemental module subroutine setMinMaxIndi_CK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if CK2_ENABLED
    pure elemental module subroutine setMinMaxIndi_CK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

#if CK1_ENABLED
    pure elemental module subroutine setMinMaxIndi_CK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(inout)             :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module subroutine setMinMaxIndi_RK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(inout)             :: a, b
    end subroutine
#endif

#if RK4_ENABLED
    pure elemental module subroutine setMinMaxIndi_RK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(inout)             :: a, b
    end subroutine
#endif

#if RK3_ENABLED
    pure elemental module subroutine setMinMaxIndi_RK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(inout)             :: a, b
    end subroutine
#endif

#if RK2_ENABLED
    pure elemental module subroutine setMinMaxIndi_RK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(inout)             :: a, b
    end subroutine
#endif

#if RK1_ENABLED
    pure elemental module subroutine setMinMaxIndi_RK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxIndi_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(inout)             :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setMinMaxPair_SK5(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(inout)             :: pair(2)
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setMinMaxPair_SK4(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(inout)             :: pair(2)
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setMinMaxPair_SK3(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(inout)             :: pair(2)
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setMinMaxPair_SK2(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(inout)             :: pair(2)
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setMinMaxPair_SK1(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(inout)             :: pair(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setMinMaxPair_IK5(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setMinMaxPair_IK4(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setMinMaxPair_IK3(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setMinMaxPair_IK2(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setMinMaxPair_IK1(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setMinMaxPair_LK5(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setMinMaxPair_LK4(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setMinMaxPair_LK3(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setMinMaxPair_LK2(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setMinMaxPair_LK1(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setMinMaxPair_CK5(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setMinMaxPair_CK4(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setMinMaxPair_CK3(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setMinMaxPair_CK2(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setMinMaxPair_CK1(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(inout)             :: pair(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setMinMaxPair_RK5(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(inout)             :: pair(2)
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setMinMaxPair_RK4(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(inout)             :: pair(2)
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setMinMaxPair_RK3(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(inout)             :: pair(2)
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setMinMaxPair_RK2(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(inout)             :: pair(2)
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setMinMaxPair_RK1(pair)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinMaxPair_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(inout)             :: pair(2)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setMinMax

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathMinMax ! LCOV_EXCL_LINE
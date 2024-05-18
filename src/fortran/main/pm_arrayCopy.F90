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
!>  This module contains procedures and generic interfaces for copying strided or indexed elements of one scalar string or vector of arbitrary
!>  intrinsic type and kind to strided or indexed elements of another scalar string or vector of the same type and kind.<br>
!>
!>  \remark
!>  The existence of this module may seem rather redundant because the Fortran standard array syntax allows easy copying array subsets.<br>
!>  However, it is impossible to perform elemental copy action for scalar strings.<br>
!>  This module provides a universal interface for copy actions for all intrinsic type and scalar strings, allowing seamless generic programming.<br>
!>
!>  \lapack{3.11}
!>  `SCOPY`, `DCOPY`, `CCOPY`, and `ZCOPY`.<br>
!>  In particular copying of subsets of scalar strings are also implemented as part of this module.
!>
!>  \see
!>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
!>  [pm_matrixInit](@ref pm_matrixInit)<br>
!>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayCopy

    use pm_kind, only: IK, LK, RK, SK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayCopy"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Copy an input scalar string or vector of arbitrary intrinsic type, kind, and size to another scalar string or vector of the same type, kind, and compatible size.<br>
    !>
    !>  \details
    !>  The functionality of this interface is readily available from the standard Fortran array syntax for all Fortran intrinsic type arrays.<br>
    !>  However, it is impossible to perform elemental copy action for scalar strings.<br>
    !>  This interface provides a universal generic approach to performing copy all vectors of intrinsic type and and kind as well as scalar strings, allowing seamless generic programming.<br>
    !>
    !>  \param[in]      From    :   The input<br>
    !>                              <ul>
    !>                                  <li>    **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ul>
    !>                              or `contiguous` array of rank `1` of either,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              whose specified elements will be copied to the array `to`.<br>
    !>  \param[inout]   To      :   The input/output scalar or `contiguous` array of the same type, kind, rank, and compatible size as the input `From`.<br>
    !>                              On output, the specific requested elements of `To` will be replaced by the requested corresponding values from `From`.<br>
    !>                              The remaining elements of `To` are returned as is.<br>
    !>  \param[in]      indexF  :   The input array of rank `1` of type `integer` of default kind \IK, of the same size as `indexT`,
    !>                              containing the indices of the elements of the input vector `From` for which copy action must be performed.<br>
    !>                              By definition, the condition `all(1_IK <= indexF) .and. all(indexF <= size(From))` must hold.<br>
    !>  \param[in]      indexT  :   The input array of rank `1` of type `integer` of default kind \IK, of the same size as `indexF`,
    !>                              containing the indices of the elements of the input/output vector `To` for which copy action must be performed.<br>
    !>                              By definition, the condition `all(1_IK <= indexT) .and. all(indexT <= size(To))` must hold.<br>
    !>
    !>  \interface{setCopyIndexed}
    !>  \code{.F90}
    !>
    !>      use pm_arrayCopy, only: setCopyIndexed
    !>
    !>      call setCopyIndexed(From, To, indexF, indexT)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `len(From) == len(To)` must hold for the corresponding scalar arguments.<br>
    !>  The condition `size(From) == size(To)` must hold for the corresponding vector arguments.<br>
    !>  The conditions `all(1_IK <= indexT) .and. all(indexT <= size(To))` must hold for the corresponding input arguments.<br>
    !>  The conditions `all(1_IK <= indexF) .and. all(indexF <= size(From))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
    !>  [pm_matrixInit](@ref pm_matrixInit)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>
    !>  \example{setCopyIndexed}
    !>  \include{lineno} example/pm_arrayCopy/setCopyIndexed/main.F90
    !>  \compilef{setCopyIndexed}
    !>  \output{setCopyIndexed}
    !>  \include{lineno} example/pm_arrayCopy/setCopyIndexed/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCopy](@ref test_pm_arrayCopy)
    !>
    !>  \final{setCopyIndexed}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface setCopyIndexed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCopyIndexed_D0_SK5(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCopyIndexed_D0_SK4(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCopyIndexed_D0_SK3(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCopyIndexed_D0_SK2(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCopyIndexed_D0_SK1(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCopyIndexed_D1_SK5(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCopyIndexed_D1_SK4(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCopyIndexed_D1_SK3(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCopyIndexed_D1_SK2(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCopyIndexed_D1_SK1(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCopyIndexed_D1_IK5(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCopyIndexed_D1_IK4(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCopyIndexed_D1_IK3(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCopyIndexed_D1_IK2(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCopyIndexed_D1_IK1(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCopyIndexed_D1_LK5(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCopyIndexed_D1_LK4(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCopyIndexed_D1_LK3(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCopyIndexed_D1_LK2(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCopyIndexed_D1_LK1(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCopyIndexed_D1_CK5(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCopyIndexed_D1_CK4(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCopyIndexed_D1_CK3(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCopyIndexed_D1_CK2(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCopyIndexed_D1_CK1(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCopyIndexed_D1_RK5(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCopyIndexed_D1_RK4(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCopyIndexed_D1_RK3(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCopyIndexed_D1_RK2(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCopyIndexed_D1_RK1(From, To, indexF, indexT)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyIndexed_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)    , contiguous    :: indexF(:), indexT(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Copy the strided elements of an input scalar string or vector of arbitrary intrinsic type, kind, and size
    !>  to the strided elements of another scalar string or vector of the same type, kind, and compatible size.<br>
    !>
    !>  \details
    !>  The functionality of this interface is readily available from the standard Fortran array syntax for all Fortran intrinsic type arrays.<br>
    !>  However, it is impossible to perform elemental copy action for scalar strings.<br>
    !>  This interface provides a universal generic approach to performing copy all vectors of intrinsic type and and kind as well as scalar strings, allowing seamless generic programming.<br>
    !>
    !>  \param[in]      From    :   The input<br>
    !>                              <ol>
    !>                                  <li>    **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              or `contiguous` array of rank `1` of either,<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter,
    !>                                  <li>    type `integer` of kind \IKALL,
    !>                                  <li>    type `logical` of kind \LKALL,
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              whose specified elements will be copied to the array `to`.<br>
    !>  \param[inout]   To      :   The input/output scalar or `contiguous` array of the same type, kind, rank, and compatible size as the input `From`.
    !>                              On output, the specific requested elements of `To` will be replaced by the requested corresponding values from `From`.
    !>                              The remaining elements of `To` are returned as is.<br>
    !>  \param[in]      incf    :   The input scalar of type `integer` of default kind \IK, containing the stride of the input vector `From`.<br>
    !>                              <ol>
    !>                                  <li>    A negative value implies the copy action to be performed from the end of `From` to its beginning every other `incf` elements.
    !>                                  <li>    A zero value will lead to copying values only from the first element of `From` to any specified element of `To`.
    !>                              </ol>
    !>  \param[in]      inct    :   The input scalar of type `integer` of default kind \IK, containing the stride of the input/output vector `To`.<br>
    !>                              The copy action will occur at every other `incf` elements in `To`.<br>
    !>                              <ol>
    !>                                  <li>    A negative value implies the copy action to be performed from the end of `To` to its beginning.
    !>                                  <li>    A zero value will lead to copying values to the first element of `To` from any specified element of `From`.
    !>                              </ol>
    !>
    !>  \interface{setCopyStrided}
    !>  \code{.F90}
    !>
    !>      use pm_arrayCopy, only: setCopyStrided
    !>
    !>      call setCopyStrided(From, To, incf, inct)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `incf /= 0_IK .or. inct /= 0_IK` must hold for the corresponding input arguments.<br>
    !>  The condition `(size(From) - 1_IK) / abs(incf) == (size(To) - 1_IK) / abs(inct)` must hold for the corresponding input arguments (when both `incf` and `inct` are non-zero).<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \lapack{3.11}
    !>  `SCOPY`, `DCOPY`, `CCOPY`, and `ZCOPY`.<br>
    !>  In particular copying of subsets of scalar strings are also implemented as part of this module.<br>
    !>
    !>  \see
    !>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
    !>  [pm_matrixInit](@ref pm_matrixInit)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>  [netlib::LAPACK](http://netlib.org/lapack/)<br>
    !>  The IBM Engineering and Scientific Subroutine Library.<br>
    !>  Developer Reference for Intel® oneAPI Math Kernel Library - Fortran.<br>
    !>  Lawson, C. L.; Hanson, R. J.; Kincaid, D. R.; Krough, F. T. Sept. 1979. “Basic Linear Algebra Subprograms for Fortran Usage.” ACM Transactions on Mathematical Software 5(3):308–323.<br>
    !>
    !>  \example{setCopyStrided}
    !>  \include{lineno} example/pm_arrayCopy/setCopyStrided/main.F90
    !>  \compilef{setCopyStrided}
    !>  \output{setCopyStrided}
    !>  \include{lineno} example/pm_arrayCopy/setCopyStrided/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCopy](@ref test_pm_arrayCopy)
    !>
    !>  \todo
    !>  \pmed
    !>  Benchmarks comparing this interface with LAPACK routines and conventional approach should be added to the documentation.<br>
    !>
    !>  \final{setCopyStrided}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface setCopyStrided

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCopyStrided_D0_SK5(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCopyStrided_D0_SK4(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCopyStrided_D0_SK3(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCopyStrided_D0_SK2(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCopyStrided_D0_SK1(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                    :: From
        character(*,SKC)        , intent(inout)                 :: To
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCopyStrided_D1_SK5(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCopyStrided_D1_SK4(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCopyStrided_D1_SK3(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCopyStrided_D1_SK2(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCopyStrided_D1_SK1(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)    , contiguous    :: From(:)
        character(*,SKC)        , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCopyStrided_D1_IK5(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCopyStrided_D1_IK4(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCopyStrided_D1_IK3(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCopyStrided_D1_IK2(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCopyStrided_D1_IK1(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , contiguous    :: From(:)
        integer(IKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCopyStrided_D1_LK5(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCopyStrided_D1_LK4(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCopyStrided_D1_LK3(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCopyStrided_D1_LK2(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCopyStrided_D1_LK1(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in)    , contiguous    :: From(:)
        logical(LKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCopyStrided_D1_CK5(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCopyStrided_D1_CK4(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCopyStrided_D1_CK3(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCopyStrided_D1_CK2(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCopyStrided_D1_CK1(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , contiguous    :: From(:)
        complex(CKC)            , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCopyStrided_D1_RK5(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCopyStrided_D1_RK4(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCopyStrided_D1_RK3(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCopyStrided_D1_RK2(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCopyStrided_D1_RK1(From, To, incf, inct)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCopyStrided_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: From(:)
        real(RKC)               , intent(inout) , contiguous    :: To(:)
        integer(IK)             , intent(in)                    :: incf, inct
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayCopy
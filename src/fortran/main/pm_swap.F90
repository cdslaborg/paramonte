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
!>  This module contains procedures and generic interfaces for swapping values of intrinsic Fortran types of arbitrary kinds.<br>
!>
!>  \details
!>  See also the procedures of [pm_arrayReverse](@ref pm_arrayReverse) for swapping the values of a vector of length `2`.<br>
!>
!>  \lapack{3.11}
!>  `SSWAP`, `DSWAP`, `CSWAP`, `ZSWAP`.<br>
!>  In particular, swapping scalars or vectors or arrays of arbitrary intrinsic 
!>  type (including `character`, `integer`, and `logical`), kind, and rank are possible.<br>
!>
!>  \see
!>  [pm_swap](@ref pm_swap)<br>
!>  [pm_arrayReverse](@ref pm_arrayReverse)<br>
!>  [pm_mathDivMul](@ref pm_mathDivMul)<br>
!>  [pm_mathSubAdd](@ref pm_mathSubAdd)<br>
!>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
!>
!>  \finmain
!>
!>  \test
!>  [test_pm_swap](@ref test_pm_swap)
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_swap

    use pm_kind, only: IK, RK, SK
    use pm_blas, only: blasSWAP

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_swap"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the two input scalar or array values `a` and `b` while their values are swapped.<br>
    !>
    !>  \param[inout]   a       :   The input/output scalar, or array of the same rank, shape, and size other array-like arguments, of either <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary `len` type-parameter, or<br>
    !>                                  <li>    type `integer` of kind \IKALL or <br>
    !>                                  <li>    type `logical` of kind \LKALL or <br>
    !>                                  <li>    type `complex` of kind \CKALL or <br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ol>
    !>                              On output, its value is swapped with the value of `b`.
    !>  \param[inout]   b       :   The input/output scalar, or array of the same rank, shape, and size other array-like arguments, of the same type and kind as the input `a`.<br>
    !>                              On output, its value is swapped with the value of `a`.
    !>  \param[in]      inca    :   The input scalar `integer` of default kind \IK, containing the stride of the input/output  scalar `character` or vector argument `a(:)`.<br>
    !>                              <ol>
    !>                                  <li>    A positive value implies the multiplication to be performed on the subset `a(1 : 1 + (size(a) - 1) * inca : inca)`.
    !>                                  <li>    A negative value implies the multiplication to be performed on the subset `a(1 + (1 - size(a)) * inca : 1 : inca)`.
    !>                              </ol>
    !>                              (**optional**, default = `1`. It can be present only if `a` and `b` are of rank `1`. It must be present if `inca` is present.)
    !>  \param[in]      incb    :   The input scalar `integer` of default kind \IK, containing the stride of the input/output scalar `character` or vector argument `b(:)`.<br>
    !>                              <ol>
    !>                                  <li>    A positive value implies the multiplication to be performed on the subset `b(1 : 1 + (size(b) - 1) * incb : incb)`.<br>
    !>                                  <li>    A negative value implies the multiplication to be performed on the subset `b(1 + (1 - size(b)) * incb : 1 : incb)`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `1`. It can be present only if `a` and `b` are of rank `1`. It must be present if `inca` is present.)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_swap, only: setSwapped
    !>
    !>      call setSwapped(a, b) ! elemental
    !>      call setSwapped(a, b, inca, incb) ! `a` and `b` are scalars of type `character`.
    !>      call setSwapped(a(..), b(..)) ! elemental
    !>      call setSwapped(a(:), b(:), inca, incb)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `len(a) == len(b)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(a(1::abs(inca))) == size(b(1::abs(incb)))` must hold for the corresponding input arguments.<br>
    !>  The rank and size of the two input array-like arguments must be the same.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  This generic interface implements a fast interface for vector input arguments.<br>
    !>  However, the cases for input arguments of higher ranks are handled via the `elemental` property of the scalar interface.<br>
    !>
    !>  \lapack{3.11}
    !>  `SSWAP`, `DSWAP`, `CSWAP`, `ZSWAP`.<br>
    !>  In particular, swapping scalars or vectors or arrays of arbitrary intrinsic 
    !>  type (including `character`, `integer`, and `logical`), kind, and rank are possible.<br>
    !>
    !>  \see
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_swap/setSwapped/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_swap/setSwapped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_swap](@ref test_pm_swap)
    !>
    !>  \todo
    !>  This generic interface can be extended to take an optional `mask` input argument that leads to selective swapping of array elements:
    !>  \verbatim
    !>  \param[in]      mask    :   The input scalar, or array of the same rank, shape, and size as other array-like arguments, of type `logical` of default kind \LK.<br>
    !>                              If present, then swapping will be done only if the corresponding value in `mask` is `.true.`.<br>
    !>                              (**optional**, default = `.true.`.)
    !>  \endverbatim
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface setSwapped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_SK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(inout)                 :: a, b
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_SK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(inout)                 :: a, b
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_SK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(inout)                 :: a, b
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_SK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(inout)                 :: a, b
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_SK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(inout)                 :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_IK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if IK4_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_IK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if IK3_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_IK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if IK2_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_IK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if IK1_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_IK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_LK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if LK4_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_LK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if LK3_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_LK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if LK2_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_LK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if LK1_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_LK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_CK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_CK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_CK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_CK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_CK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(inout)                 :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_RK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(inout)                 :: a, b
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_RK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(inout)                 :: a, b
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_RK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(inout)                 :: a, b
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_RK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(inout)                 :: a, b
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setSwappedDef_D0_RK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(inout)                 :: a, b
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSwappedDef_D1_SK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSwappedDef_D1_SK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSwappedDef_D1_SK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSwappedDef_D1_SK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSwappedDef_D1_SK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSwappedDef_D1_IK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSwappedDef_D1_IK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSwappedDef_D1_IK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSwappedDef_D1_IK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSwappedDef_D1_IK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSwappedDef_D1_LK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSwappedDef_D1_LK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSwappedDef_D1_LK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSwappedDef_D1_LK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSwappedDef_D1_LK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSwappedDef_D1_CK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSwappedDef_D1_CK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSwappedDef_D1_CK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSwappedDef_D1_CK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSwappedDef_D1_CK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSwappedDef_D1_RK5(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSwappedDef_D1_RK4(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSwappedDef_D1_RK3(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSwappedDef_D1_RK2(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSwappedDef_D1_RK1(a, b)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedDef_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSwappedInc_D0_SK5(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(inout)                 :: a, b
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSwappedInc_D0_SK4(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(inout)                 :: a, b
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSwappedInc_D0_SK3(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(inout)                 :: a, b
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSwappedInc_D0_SK2(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(inout)                 :: a, b
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSwappedInc_D0_SK1(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(inout)                 :: a, b
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSwappedInc_D1_SK5(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSwappedInc_D1_SK4(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSwappedInc_D1_SK3(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSwappedInc_D1_SK2(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSwappedInc_D1_SK1(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                        , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSwappedInc_D1_IK5(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSwappedInc_D1_IK4(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSwappedInc_D1_IK3(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSwappedInc_D1_IK2(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSwappedInc_D1_IK1(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSwappedInc_D1_LK5(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSwappedInc_D1_LK4(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSwappedInc_D1_LK3(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSwappedInc_D1_LK2(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSwappedInc_D1_LK1(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSwappedInc_D1_CK5(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSwappedInc_D1_CK4(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSwappedInc_D1_CK3(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSwappedInc_D1_CK2(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSwappedInc_D1_CK1(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                            , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSwappedInc_D1_RK5(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSwappedInc_D1_RK4(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSwappedInc_D1_RK3(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSwappedInc_D1_RK2(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSwappedInc_D1_RK1(a, b, inca, incb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSwappedInc_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                               , intent(inout) , contiguous    :: a(:), b(:)
        integer(IK)                             , intent(in)                    :: inca, incb
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setSwapped

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_swap ! LCOV_EXCL_LINE
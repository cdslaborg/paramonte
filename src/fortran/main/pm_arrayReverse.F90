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
!>  This module contains procedures and generic interfaces for reversing the order of elements in arrays of various types.
!>
!>  \benchmarks
!>
!>  \benchmark{getReversed_vs_setReversed, The runtime performance of [setReversed](@ref pm_arrayReverse::setReversed) vs. direct reversing}
!>  \include{lineno} benchmark/pm_arrayReverse/getReversed_vs_setReversed/main.F90
!>  \compilefb{getReversed_vs_setReversed}
!>  \postprocb{getReversed_vs_setReversed}
!>  \include{lineno} benchmark/pm_arrayReverse/getReversed_vs_setReversed/main.py
!>  \visb{getReversed_vs_setReversed}
!>  \image html benchmark/pm_arrayReverse/getReversed_vs_setReversed/benchmark.getReversed_vs_setReversed.runtime.png width=1000
!>  \image html benchmark/pm_arrayReverse/getReversed_vs_setReversed/benchmark.getReversed_vs_setReversed.runtime.ratio.png width=1000
!>  \moralb{getReversed_vs_setReversed}
!>      -#  The procedures under the generic interface [setReversed](@ref pm_arrayReverse::setReversed) with in-place array reversal
!>          (corresponding to `setReversed_overwrite()` in the benchmark) tend to be significantly faster than the functional-interface
!>          procedures [getReversed](@ref pm_arrayReverse::getReversed).<br>
!>          The reason is likely the extra copy required with the functional interface.
!>
!>  \test
!>  [test_pm_arrayReverse](@ref test_pm_arrayReverse)
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
module pm_arrayReverse

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayReverse"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an output array whose elements are the reversed-order elements of the input `Array`, such that <br>
    !>  &nbsp; `Array = Array(lenArray:1:-1)`, <br>
    !>  &nbsp; where `lenArray = len(Array)` if `array` is a scalar character, or `lenArray = size(Array)` is an array of rank `1`.
    !>
    !>  \param[in]  Array   :   The input `contiguous` array of shape `(:)` of either <br>
    !>                          <ul>
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                              <li>    type `character` of kind \SKALL, or <br>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                              <li>    type `logical` of kind \LKALL, or <br>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a scalar `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                          </ul>
    !>                          whose elements will be reversed in order.<br>
    !>
    !>  \return
    !>  `ArrayReversed`     :   The output `contiguous` array of the same type, kind, shape, and size as the input `array` that will contain the reversed array.
    !>
    !>  \interface{getReversed}
    !>  \code{.F90}
    !>
    !>      use pm_arrayReverse, only: getReversed
    !>
    !>      ArrayReversed = getReversed(Array)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The primary purpose of the procedures under this generic interface is to provide an efficient **generic** method of reversing an array.<br>
    !>  This generic interface is particularly useful in case of scalar character input arguments.<br>
    !>  With further compiler and language template enhancements in the future,
    !>  the need for the procedures under this generic interface might be resolved in the future.<br>
    !>  See [pm_arrayReverse](@ref pm_arrayReverse) for the relevant benchmarks.<br>
    !>
    !>  \see
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>
    !>  \example{getReversed}
    !>  \include{lineno} example/pm_arrayReverse/getReversed/main.F90
    !>  \compilef{getReversed}
    !>  \output{getReversed}
    !>  \include{lineno} example/pm_arrayReverse/getReversed/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayReverse](@ref test_pm_arrayReverse)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \finmain{getReversed}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface getReversed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReversedNew_D0_SK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                    :: Array
        character(len(Array,IK),SKC)                                :: ArrayReversed
    end function
#endif

#if SK4_ENABLED
    PURE module function getReversedNew_D0_SK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                    :: Array
        character(len(Array,IK),SKC)                                :: ArrayReversed
    end function
#endif

#if SK3_ENABLED
    PURE module function getReversedNew_D0_SK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                    :: Array
        character(len(Array,IK),SKC)                                :: ArrayReversed
    end function
#endif

#if SK2_ENABLED
    PURE module function getReversedNew_D0_SK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                    :: Array
        character(len(Array,IK),SKC)                                :: ArrayReversed
    end function
#endif

#if SK1_ENABLED
    PURE module function getReversedNew_D0_SK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                    :: Array
        character(len(Array,IK),SKC)                                :: ArrayReversed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReversedNew_D1_SK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(len(Array,IK),SKC)                                :: ArrayReversed(size(Array))
    end function
#endif

#if SK4_ENABLED
    PURE module function getReversedNew_D1_SK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(len(Array,IK),SKC)                                :: ArrayReversed(size(Array))
    end function
#endif

#if SK3_ENABLED
    PURE module function getReversedNew_D1_SK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(len(Array,IK),SKC)                                :: ArrayReversed(size(Array))
    end function
#endif

#if SK2_ENABLED
    PURE module function getReversedNew_D1_SK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(len(Array,IK),SKC)                                :: ArrayReversed(size(Array))
    end function
#endif

#if SK1_ENABLED
    PURE module function getReversedNew_D1_SK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(len(Array,IK),SKC)                                :: ArrayReversed(size(Array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReversedNew_D1_IK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if IK4_ENABLED
    PURE module function getReversedNew_D1_IK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if IK3_ENABLED
    PURE module function getReversedNew_D1_IK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if IK2_ENABLED
    PURE module function getReversedNew_D1_IK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if IK1_ENABLED
    PURE module function getReversedNew_D1_IK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReversedNew_D1_LK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if LK4_ENABLED
    PURE module function getReversedNew_D1_LK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if LK3_ENABLED
    PURE module function getReversedNew_D1_LK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if LK2_ENABLED
    PURE module function getReversedNew_D1_LK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if LK1_ENABLED
    PURE module function getReversedNew_D1_LK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReversedNew_D1_CK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if CK4_ENABLED
    PURE module function getReversedNew_D1_CK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if CK3_ENABLED
    PURE module function getReversedNew_D1_CK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if CK2_ENABLED
    PURE module function getReversedNew_D1_CK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

#if CK1_ENABLED
    PURE module function getReversedNew_D1_CK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                                                :: ArrayReversed(size(Array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReversedNew_D1_RK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                                                   :: ArrayReversed(size(Array))
    end function
#endif

#if RK4_ENABLED
    PURE module function getReversedNew_D1_RK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                                                   :: ArrayReversed(size(Array))
    end function
#endif

#if RK3_ENABLED
    PURE module function getReversedNew_D1_RK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                                                   :: ArrayReversed(size(Array))
    end function
#endif

#if RK2_ENABLED
    PURE module function getReversedNew_D1_RK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                                                   :: ArrayReversed(size(Array))
    end function
#endif

#if RK1_ENABLED
    PURE module function getReversedNew_D1_RK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                                                   :: ArrayReversed(size(Array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReversedNew_D1_PSSK5(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_PSSK5
#endif
        use pm_kind, only: SKC => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))                                       :: ArrayReversed(size(Array))
    end function
#endif

#if SK4_ENABLED
    PURE module function getReversedNew_D1_PSSK4(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_PSSK4
#endif
        use pm_kind, only: SKC => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))                                       :: ArrayReversed(size(Array))
    end function
#endif

#if SK3_ENABLED
    PURE module function getReversedNew_D1_PSSK3(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_PSSK3
#endif
        use pm_kind, only: SKC => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))                                       :: ArrayReversed(size(Array))
    end function
#endif

#if SK2_ENABLED
    PURE module function getReversedNew_D1_PSSK2(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_PSSK2
#endif
        use pm_kind, only: SKC => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))                                       :: ArrayReversed(size(Array))
    end function
#endif

#if SK1_ENABLED
    PURE module function getReversedNew_D1_PSSK1(Array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversed_D1_PSSK1
#endif
        use pm_kind, only: SKC => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))                                       :: ArrayReversed(size(Array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Reverse the order of the elements of the input `Array`, such that <br>
    !>  &nbsp; `Array = Array(lenArray:1:-1)`, <br>
    !>  &nbsp; where `lenArray = len(Array)` if `array` is a scalar character, or `lenArray = size(Array)` is an array of rank `1`.
    !>
    !>  \param[inout]   Array       :   The **input** or **input/output** `contiguous` array of shape `(:)` of either <br>
    !>                                  <ul>
    !>                                      <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                                      <li>    type `character` of kind \SKALL, or <br>
    !>                                      <li>    type `integer` of kind \IKALL, or <br>
    !>                                      <li>    type `logical` of kind \LKALL, or <br>
    !>                                      <li>    type `complex` of kind \CKALL, or <br>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ul>
    !>                                  or,
    !>                                  <ul>
    !>                                      <li>    a scalar `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                                  </ul>
    !>                                  whose elements will be reversed in order, such that the first element becomes the last and the last becomes the first in Array.<br>
    !>                                  <ul>
    !>                                      <li>    If the input argument `ArrayReversed` is specified, then `array` has `intent(in)` and will not change.<br>
    !>                                      <li>    If the input argument `ArrayReversed` is missing, then `array` will be reversed on output.<br>
    !>                                  </ul>
    !>  \param[out] ArrayReversed   :   The output `contiguous` array of the same type, kind, shape, and size as the input `array` that will contain the reversed array.<br>
    !>                                  (**optional**, if missing, the reversed array will stored and output in the input contiguous `array` **in-place**)
    !>
    !>  \interface{setReversed}
    !>  \code{.F90}
    !>
    !>      use pm_arrayReverse, only: setReversed
    !>
    !>      call setReversed(Array)
    !>      call setReversed(Array, ArrayReversed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If `ArrayReversed` is present, its size must equal that of `Array`.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The primary purpose of the procedures under this generic interface is to provide a more efficient faster method of reversing an,
    !>  array of any intrinsic type and kind, in place, without an extra copy that is implicitly done by current compilers with the regular Fortran syntax.<br>
    !>  This generic interface also provides a comfortable generic way to reverse arrays of any intrinsic type and kind. This is particularly useful in case of scalar character.<br>
    !>  With further compiler and language template enhancements in the future, the need for the procedures under this generic interface might be resolved in the future.<br>
    !>  See [pm_arrayReverse](@ref pm_arrayReverse) for the relevant benchmarks.<br>
    !>
    !>  \note
    !>  Upon return, if `array` is `allocatable, intent(inout) :: Array` argument, then it is guaranteed to have the same lower and upper bounds as before.
    !>  This condition happens when `ArrayReversed` is missing in the arguments list.
    !>
    !>  \see
    !>  [getReversed](@ref pm_arrayReverse::getReversed)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>
    !>  \example{setReversed}
    !>  \include{lineno} example/pm_arrayReverse/setReversed/main.F90
    !>  \compilef{setReversed}
    !>  \output{setReversed}
    !>  \include{lineno} example/pm_arrayReverse/setReversed/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayReverse](@ref test_pm_arrayReverse)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \finmain{setReversed}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface setReversed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedOld_D0_SK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(inout)                 :: Array
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedOld_D0_SK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(inout)                 :: Array
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedOld_D0_SK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(inout)                 :: Array
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedOld_D0_SK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(inout)                 :: Array
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedOld_D0_SK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(inout)                 :: Array
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedOld_D1_SK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedOld_D1_SK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedOld_D1_SK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedOld_D1_SK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedOld_D1_SK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReversedOld_D1_IK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReversedOld_D1_IK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReversedOld_D1_IK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReversedOld_D1_IK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setReversedOld_D1_IK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReversedOld_D1_LK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReversedOld_D1_LK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReversedOld_D1_LK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReversedOld_D1_LK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setReversedOld_D1_LK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReversedOld_D1_CK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReversedOld_D1_CK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReversedOld_D1_CK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReversedOld_D1_CK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setReversedOld_D1_CK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReversedOld_D1_RK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReversedOld_D1_RK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReversedOld_D1_RK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReversedOld_D1_RK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setReversedOld_D1_RK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK5(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_PSSK5
#endif
        use pm_kind, only: SKC => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK4(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_PSSK4
#endif
        use pm_kind, only: SKC => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK3(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_PSSK3
#endif
        use pm_kind, only: SKC => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK2(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_PSSK2
#endif
        use pm_kind, only: SKC => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK1(Array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversed_D1_PSSK1
#endif
        use pm_kind, only: SKC => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(inout) , contiguous    :: Array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedNew_D0_SK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                    :: Array
        character(*,SKC)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedNew_D0_SK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                    :: Array
        character(*,SKC)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedNew_D0_SK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                    :: Array
        character(*,SKC)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedNew_D0_SK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                    :: Array
        character(*,SKC)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedNew_D0_SK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                    :: Array
        character(*,SKC)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedNew_D1_SK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(*,SKC)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedNew_D1_SK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(*,SKC)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedNew_D1_SK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(*,SKC)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedNew_D1_SK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(*,SKC)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedNew_D1_SK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: Array(:)
        character(*,SKC)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReversedNew_D1_IK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReversedNew_D1_IK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReversedNew_D1_IK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReversedNew_D1_IK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setReversedNew_D1_IK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: Array(:)
        integer(IKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReversedNew_D1_LK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReversedNew_D1_LK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReversedNew_D1_LK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReversedNew_D1_LK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setReversedNew_D1_LK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: Array(:)
        logical(LKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReversedNew_D1_CK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReversedNew_D1_CK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReversedNew_D1_CK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReversedNew_D1_CK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setReversedNew_D1_CK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: Array(:)
        complex(CKC)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReversedNew_D1_RK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReversedNew_D1_RK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReversedNew_D1_RK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReversedNew_D1_RK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setReversedNew_D1_RK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: Array(:)
        real(RKC)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK5(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK5
#endif
        use pm_kind, only: SKC => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK4(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK4
#endif
        use pm_kind, only: SKC => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK3(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK3
#endif
        use pm_kind, only: SKC => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK2(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK2
#endif
        use pm_kind, only: SKC => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK1(Array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK1
#endif
        use pm_kind, only: SKC => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKC))       , intent(in)    , contiguous    :: Array(:)
        type(css_pdt(SKC))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayReverse ! LCOV_EXCL_LINE
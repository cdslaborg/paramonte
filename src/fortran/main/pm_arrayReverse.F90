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
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
module pm_arrayReverse

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayReverse"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an output array whose elements are the reversed-order elements of the input `array`, such that <br>
    !>  &nbsp; `array = array(lenArray:1:-1)`, <br>
    !>  &nbsp; where `lenArray = len(array)` if `array` is a scalar character, or `lenArray = size(array)` is an array of rank `1`.
    !>
    !>  \param[in]  array   :   The input `contiguous` array of shape `(:)` of either <br>
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
    !>      ArrayReversed = getReversed(array)
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
    !>  \final{getReversed}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getReversed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReversedNew_D0_SK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(len(array,IK),SKG)                                :: ArrayReversed
    end function
#endif

#if SK4_ENABLED
    PURE module function getReversedNew_D0_SK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(len(array,IK),SKG)                                :: ArrayReversed
    end function
#endif

#if SK3_ENABLED
    PURE module function getReversedNew_D0_SK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(len(array,IK),SKG)                                :: ArrayReversed
    end function
#endif

#if SK2_ENABLED
    PURE module function getReversedNew_D0_SK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(len(array,IK),SKG)                                :: ArrayReversed
    end function
#endif

#if SK1_ENABLED
    PURE module function getReversedNew_D0_SK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(len(array,IK),SKG)                                :: ArrayReversed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReversedNew_D1_SK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                                :: ArrayReversed(size(array))
    end function
#endif

#if SK4_ENABLED
    PURE module function getReversedNew_D1_SK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                                :: ArrayReversed(size(array))
    end function
#endif

#if SK3_ENABLED
    PURE module function getReversedNew_D1_SK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                                :: ArrayReversed(size(array))
    end function
#endif

#if SK2_ENABLED
    PURE module function getReversedNew_D1_SK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                                :: ArrayReversed(size(array))
    end function
#endif

#if SK1_ENABLED
    PURE module function getReversedNew_D1_SK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                                :: ArrayReversed(size(array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReversedNew_D1_IK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if IK4_ENABLED
    PURE module function getReversedNew_D1_IK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if IK3_ENABLED
    PURE module function getReversedNew_D1_IK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if IK2_ENABLED
    PURE module function getReversedNew_D1_IK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if IK1_ENABLED
    PURE module function getReversedNew_D1_IK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                                :: ArrayReversed(size(array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReversedNew_D1_LK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if LK4_ENABLED
    PURE module function getReversedNew_D1_LK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if LK3_ENABLED
    PURE module function getReversedNew_D1_LK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if LK2_ENABLED
    PURE module function getReversedNew_D1_LK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if LK1_ENABLED
    PURE module function getReversedNew_D1_LK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                                :: ArrayReversed(size(array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReversedNew_D1_CK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if CK4_ENABLED
    PURE module function getReversedNew_D1_CK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if CK3_ENABLED
    PURE module function getReversedNew_D1_CK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if CK2_ENABLED
    PURE module function getReversedNew_D1_CK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                                :: ArrayReversed(size(array))
    end function
#endif

#if CK1_ENABLED
    PURE module function getReversedNew_D1_CK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                                :: ArrayReversed(size(array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReversedNew_D1_RK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                                   :: ArrayReversed(size(array))
    end function
#endif

#if RK4_ENABLED
    PURE module function getReversedNew_D1_RK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                                   :: ArrayReversed(size(array))
    end function
#endif

#if RK3_ENABLED
    PURE module function getReversedNew_D1_RK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                                   :: ArrayReversed(size(array))
    end function
#endif

#if RK2_ENABLED
    PURE module function getReversedNew_D1_RK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                                   :: ArrayReversed(size(array))
    end function
#endif

#if RK1_ENABLED
    PURE module function getReversedNew_D1_RK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                                   :: ArrayReversed(size(array))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getReversedNew_D1_PSSK5(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                       :: ArrayReversed(size(array))
    end function
#endif

#if SK4_ENABLED
    PURE module function getReversedNew_D1_PSSK4(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                       :: ArrayReversed(size(array))
    end function
#endif

#if SK3_ENABLED
    PURE module function getReversedNew_D1_PSSK3(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                       :: ArrayReversed(size(array))
    end function
#endif

#if SK2_ENABLED
    PURE module function getReversedNew_D1_PSSK2(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                       :: ArrayReversed(size(array))
    end function
#endif

#if SK1_ENABLED
    PURE module function getReversedNew_D1_PSSK1(array) result(ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReversedNew_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                       :: ArrayReversed(size(array))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Reverse the order of the elements of the input `array`, such that <br>
    !>  &nbsp; `array = array(lenArray:1:-1)`, <br>
    !>  &nbsp; where `lenArray = len(array)` if `array` is a scalar character, or `lenArray = size(array)` is an array of rank `1`.
    !>
    !>  \param[inout]   array       :   The **input** or **input/output** `contiguous` array of shape `(:)` of either <br>
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
    !>                                  whose elements will be reversed in order, such that the first element becomes the last and the last becomes the first in array.<br>
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
    !>      call setReversed(array)
    !>      call setReversed(array, ArrayReversed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If `ArrayReversed` is present, its size must equal that of `array`.<br>
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
    !>  Upon return, if `array` is `allocatable, intent(inout) :: array` argument, then it is guaranteed to have the same lower and upper bounds as before.
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
    !>  \final{setReversed}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setReversed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedOld_D0_SK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedOld_D0_SK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedOld_D0_SK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedOld_D0_SK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedOld_D0_SK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedOld_D1_SK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedOld_D1_SK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedOld_D1_SK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedOld_D1_SK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedOld_D1_SK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReversedOld_D1_IK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReversedOld_D1_IK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReversedOld_D1_IK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReversedOld_D1_IK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setReversedOld_D1_IK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReversedOld_D1_LK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReversedOld_D1_LK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReversedOld_D1_LK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReversedOld_D1_LK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setReversedOld_D1_LK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReversedOld_D1_CK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReversedOld_D1_CK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReversedOld_D1_CK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReversedOld_D1_CK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setReversedOld_D1_CK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReversedOld_D1_RK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReversedOld_D1_RK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReversedOld_D1_RK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReversedOld_D1_RK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setReversedOld_D1_RK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedOld_D1_PSSK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedOld_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedNew_D0_SK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedNew_D0_SK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedNew_D0_SK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedNew_D0_SK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedNew_D0_SK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(out)                   :: ArrayReversed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReversedNew_D1_SK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedNew_D1_SK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedNew_D1_SK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedNew_D1_SK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedNew_D1_SK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReversedNew_D1_IK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReversedNew_D1_IK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReversedNew_D1_IK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReversedNew_D1_IK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setReversedNew_D1_IK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReversedNew_D1_LK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReversedNew_D1_LK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReversedNew_D1_LK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReversedNew_D1_LK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setReversedNew_D1_LK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReversedNew_D1_CK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReversedNew_D1_CK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReversedNew_D1_CK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReversedNew_D1_CK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setReversedNew_D1_CK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReversedNew_D1_RK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReversedNew_D1_RK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReversedNew_D1_RK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReversedNew_D1_RK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setReversedNew_D1_RK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK5(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK4(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK3(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK2(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setReversedNew_D1_PSSK1(array, ArrayReversed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReversedNew_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))       , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))       , intent(out)   , contiguous    :: ArrayReversed(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayReverse ! LCOV_EXCL_LINE